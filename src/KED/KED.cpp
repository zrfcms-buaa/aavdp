#include "KED.h"

void copy_knode_data(KED_KNODE *knode1, KED_KNODE *knode2)
{
    vector_copy(knode1->hkl, knode2->hkl);
    vector_copy(knode1->K, knode2->K);
    knode1->Kmagnitude=knode2->Kmagnitude;
    knode1->intensity=knode2->intensity;
}

void swap_knode_data(KED_KNODE *knode1, KED_KNODE *knode2)
{
    KED_KNODE *ktemp=new KED_KNODE;
    copy_knode_data(ktemp, knode1);
    copy_knode_data(knode1, knode2);
    copy_knode_data(knode2, ktemp);
    delete ktemp;
}

void quick_sort(KED_KNODE *kstart, KED_KNODE *kend)
{
    if(kstart==nullptr||kend==nullptr||kstart==kend) return;
    KED_KNODE *knode1=kstart;
    KED_KNODE *knode2=kstart->next;
    double Kmagnitude=kstart->Kmagnitude;
    while(knode2!=kend->next&&knode2!=nullptr){
        if(knode2->Kmagnitude<Kmagnitude){
            knode1=knode1->next;
            if(knode1!=knode2){
                swap_knode_data(knode1, knode2);
            }
        }
        knode2=knode2->next;
    }
    swap_knode_data(knode1, kstart);
    quick_sort(kstart, knode1);
    quick_sort(knode1->next, kend);
}

void quick_sort(KED_KNODE *kstart)
{
    if(kstart==nullptr) return;
    KED_KNODE *knode1=kstart;
    KED_KNODE *knode2=kstart->next;
    double hkl[3]; vector_copy(hkl, kstart->hkl);
    while(knode2!=nullptr){
        if((0==int(knode2->hkl[0])+int(hkl[0]))&&(0==int(knode2->hkl[1])+int(hkl[1]))&&(0==int(knode2->hkl[2])+int(hkl[2]))){
            knode1=knode1->next;
            if(knode1!=knode2){
                swap_knode_data(knode1, knode2);
            }
            break;
        }else{
            knode2=knode2->next;
        }
        
    }
    quick_sort(knode1->next);
}

KED::KED(MODEL *model, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if(mpi_rank==0){
        printf("[INFO] Starting computation of kinematic electron diffraction...\n");
        printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", model->lambda);
    }
    lambda=model->lambda;
    Kmagnitude_max=Kmag_max;
    double spacingK[3];
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    int kmin[3], kmax[3];
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    int num=(2*kmax[0]+1)*(2*kmax[1]+1)*(2*kmax[2]+1);

    clock_t start, finish;
    start=clock();
    int my_count=0, task_id=0;
    if(mpi_rank==0){
        printf("[INFO] Spacings along a*, b*, and c* in reciprocal space [Angstrom-1]: %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
        printf("[INFO] Number of spacings along a*, b*, and c* in reciprocal space: %d %d %d\n", kmax[0], kmax[1], kmax[2]);
        printf("[INFO] Starting computation of diffraction intensity with %d k-points and %d processes ...\n", num, mpi_size);
        double hkl0[3]={0.0}, K0[3]={0.0};
        double intensity0=model->get_diffraction_intensity(0.0, K0, true);
        add_k_node(hkl0, K0, 0.0, intensity0); my_count++;
    }
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                task_id++;
                if(task_id%mpi_size!=mpi_rank) continue;
                if(0==ih&&0==ik&&0==il) continue;
                double hkl[3]={double(ih), double(ik), double(il)};
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmagnitude=model->get_reciprocal_vector_length(K);
                if(Kmagnitude<Kmagnitude_max){
                    model->reciprocal_to_cartesian(K, K);
                    double intensity=model->get_diffraction_intensity(Kmagnitude, K, false);
                    if(intensity>KED_INTENSITY_LIMIT){
                        add_k_node(hkl, K, Kmagnitude, intensity);
                    }
                }
                my_count++;
            }
        }
    }
    finish=clock();
    double my_time=double(finish-start)/CLOCKS_PER_SEC;
    int total_count=0;
    MPI_Reduce(&my_count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double total_time=0.0;
    MPI_Reduce(&my_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    merge_k_node();
    if(mpi_rank==0){
        printf("[INFO] Ending computation of diffraction intensity with %d k-points and %d processes\n", total_count, mpi_size);
        printf("[INFO] Computation time [s]: %.8f\n", total_time);
        printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
        printf("[INFO] Total number of diffraction intensity (including intensity at the transmission spot): %d\n", numk);
        printf("[INFO] Total range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
        filter_diffraction_intensity(threshold);
        quick_sort(khead, ktail);
        printf("[INFO] Total number of filtered diffraction intensity (including intensity at the transmission spot): %d\n", numk);
        printf("[INFO] Total range of filtered diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
        printf("[INFO] Ending computation of kinematic electron diffraction\n");
    }
}

KED::KED(MODEL *model, int zone[3], double thickness, double Kmag_max, double threshold, double spacing[3], bool is_spacing_auto)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if(mpi_rank==0){
        printf("[INFO] Starting computation of kinematic electron diffraction...\n");
        printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", model->lambda);
    }
    lambda=model->lambda;
    Kmagnitude_max=Kmag_max;
    double spacingK[3];
    if(is_spacing_auto){
        model->compute_reciprocal_spacing(spacingK, spacing);
    }else{
        vector_copy(spacingK, spacing);
    }
    int NspacingK[3];
    for(int i=0;i<3;i++){
        NspacingK[i]=ceil(Kmagnitude_max/spacingK[i]);
    }
    int kmin[3], kmax[3];
    vector_copy(kmax, NspacingK); vector_constant(kmin, -1, NspacingK);
    double n_zone[3]={double(zone[0]), double(zone[1]), double(zone[2])};
    vector_normalize(n_zone, n_zone);
    double upper_bound=thickness/2.0;
    double lower_bound=-thickness/2.0;
    int num=(2*kmax[0]+1)*(2*kmax[1]+1)*(2*kmax[2]+1);

    clock_t start, finish;
    start=clock();
    int my_count=0, task_id=0;
    if(mpi_rank==0){
        printf("[INFO] Spacings along a*, b*, and c* in reciprocal space [Angstrom-1]: %.8f %.8f %.8f\n", spacingK[0], spacingK[1], spacingK[2]);
        printf("[INFO] Number of spacings along a*, b*, and c* in reciprocal space: %d %d %d\n", kmax[0], kmax[1], kmax[2]);
        printf("[INFO] Range along zone-[%.8f %.8f %.8f] in reciprocal space [Angstrom-1]: %.8f %.8f\n", n_zone[0], n_zone[1], n_zone[2], lower_bound, upper_bound);
        printf("[INFO] Starting computation of diffraction intensity with %d k-points and %d processes ...\n", num, mpi_size);
        double hkl0[3]={0.0}, K0[3]={0.0};
        double intensity0=model->get_diffraction_intensity(0.0, K0, true);
        add_k_node(hkl0, K0, 0.0, intensity0); my_count++;
    }
    for(int ih=kmin[0];ih<=kmax[0];ih++){
        for(int ik=kmin[1];ik<=kmax[1];ik++){
            for(int il=kmin[2];il<=kmax[2];il++){
                task_id++;
                if(task_id%mpi_size!=mpi_rank) continue;
                if(0==ih&&0==ik&&0==il) continue;
                double hkl[3]={double(ih), double(ik), double(il)};
                double K[3]={double(ih)*spacingK[0], double(ik)*spacingK[1], double(il)*spacingK[2]};
                double Kmagnitude=model->get_reciprocal_vector_length(K);
                if(Kmagnitude<Kmagnitude_max){
                    model->reciprocal_to_cartesian(K, K);
                    double proj=vector_dot(K, n_zone);
                    if((proj>lower_bound)&&(proj<upper_bound)){
                        double intensity=model->get_diffraction_intensity(Kmagnitude, K, false);
                        if(intensity>KED_INTENSITY_LIMIT){
                            add_k_node(hkl, K, Kmagnitude, intensity);
                        }
                    }
                }
                my_count++;
            }
        }
    }
    finish=clock();
    double my_time=double(finish-start)/CLOCKS_PER_SEC;
    int total_count=0;
    MPI_Reduce(&my_count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double total_time=0.0;
    MPI_Reduce(&my_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    merge_k_node();
    if(mpi_rank==0){
        printf("[INFO] Ending computation of diffraction intensity with %d k-points and %d processes\n", total_count, mpi_size);
        printf("[INFO] Computation time [s]: %.8f\n", total_time);
        printf("[INFO] Intensity at the transmission spot: %.8f\n", khead->intensity);
        printf("[INFO] Total number  of diffraction intensity (including intensity at the transmission spot): %d\n", numk);
        printf("[INFO] Total range of diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
        filter_diffraction_intensity(threshold);
        printf("[INFO] Total number of filtered diffraction intensity (including intensity at the transmission spot): %d\n", numk);
        printf("[INFO] Total range of filtered diffraction intensity: %.8f %.8f\n", intensity_min, intensity_max);
        find_first_and_second_knearests();
        if(knearest_2==nullptr){
            printf("[INFO] The first nearest diffraction vectors R1: [%.8f %.8f %.8f]\n", knearest_1->K[0], knearest_1->K[1], knearest_1->K[2]);
            printf("[WARN] Unable to find the second nearest diffraction vector\n");
        }else{
            printf("[INFO] The first and second nearest diffraction vectors R1, R2: [%.8f %.8f %.8f], [%.8f %.8f %.8f] (R2/R1 %.8f and angle %.8f)\n", 
            knearest_1->K[0], knearest_1->K[1], knearest_1->K[2], knearest_2->K[0], knearest_2->K[1], knearest_2->K[2], knearest_2->Kmagnitude/knearest_1->Kmagnitude, vector_angle(knearest_2->K, knearest_1->K)*RAD_TO_DEG);
        }
        rotate_by_first_knearest(zone);
        printf("[INFO] Ending computation of kinematic electron diffraction\n");
    }
}

KED::~KED()
{
    free_k_node();
}

void KED::add_k_node(double hkl[3], double K[3], double Kmagnitude, double intensity)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new KED_KNODE;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kmagnitude=Kmagnitude;
        ktail->intensity=intensity;
        numk++;
    }else{
        ktail->next=new KED_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kmagnitude=Kmagnitude;
        ktail->intensity=intensity;
        if(intensity_max<intensity) intensity_max=intensity;
        if(intensity_min>intensity) intensity_min=intensity;
        numk++;
    }
}

void KED::merge_k_node()
{
    struct NODE{double h,k,l; double K1, K2, K3; double Kmagnitude; double intensity;};
    MPI_Datatype MPI_NODE;
    int blocklengths[4] = {3, 3, 1, 1};
    MPI_Aint offsets[4] = {offsetof(NODE, h), offsetof(NODE, K1), offsetof(NODE, Kmagnitude), offsetof(NODE, intensity)};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(4, blocklengths, offsets, types, &MPI_NODE);
    MPI_Type_commit(&MPI_NODE);

    int numk_all=0;
    MPI_Reduce(&numk, &numk_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    std::vector<int> numks(mpi_size);
    MPI_Gather(&numk,1,MPI_INT,numks.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    std::vector<NODE> vec;
    KED_KNODE *ktemp=khead;
    while(ktemp){
        vec.push_back({ktemp->hkl[0],ktemp->hkl[1],ktemp->hkl[2],ktemp->K[0],ktemp->K[1],ktemp->K[2],ktemp->Kmagnitude,ktemp->intensity});
        ktemp=ktemp->next;
    }

    std::vector<int> disps(mpi_size);
    std::vector<NODE> vec_all;
    if(mpi_rank==0){
        vec_all.resize(numk_all);
        disps[0]=0;
        for(int i=1; i<mpi_size; ++i){
            disps[i]=disps[i-1]+numks[i-1];
        }
    }
    MPI_Gatherv(vec.data(), numk, MPI_NODE, vec_all.data(), numks.data(), disps.data(), MPI_NODE, 0, MPI_COMM_WORLD);
    free_k_node();
    
    if(mpi_rank==0){
        for(auto &vec : vec_all){
            double hkl[3]={vec.h,vec.k,vec.l};
            double K[3]={vec.K1,vec.K2,vec.K3};
            add_k_node(hkl, K, vec.Kmagnitude, vec.intensity);
        }
    }
    MPI_Type_free(&MPI_NODE);
}

void KED::free_k_node()
{
    KED_KNODE *cur=khead;
    while(cur!=nullptr){
        KED_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
    khead=ktail=nullptr;
    numk=0;  
}

void KED::filter_diffraction_intensity(double threshold)
{
    double intensity_threshold=intensity_max*threshold;
    intensity_min=1.0e8;
    while(khead!=nullptr&&khead->intensity<=intensity_threshold){
        KED_KNODE *ktemp=khead;
        khead=khead->next;
        delete ktemp;
        numk--;
    }
    KED_KNODE *cur=khead;
    while(cur!=nullptr&&cur->next!=nullptr){
        if(cur->next->intensity<=intensity_threshold){
            KED_KNODE *ktemp=cur->next;
            cur->next=cur->next->next;
            delete ktemp;
            numk--;
        }else{
            if(intensity_min>cur->intensity) intensity_min=cur->intensity;
            cur=cur->next;
        }
    }
}

void KED::find_first_and_second_knearests()
{
    quick_sort(khead, ktail);
    knearest_1=khead->next;
    if(knearest_1==nullptr){
        printf("[ERROR] Unable to find the first nearest diffraction vector\n");
        exit(1);
    }
    knearest_2=knearest_1->next;
    while(knearest_2!=nullptr){
        if(knearest_2->Kmagnitude-knearest_1->Kmagnitude>KED_KMAG_LIMIT&&vector_angle(knearest_2->K, knearest_1->K)<=HALF_PI&&vector_angle(knearest_2->K, knearest_1->K)>1.0e-6) break;
        knearest_2=knearest_2->next;
    }
}

void KED::rotate_by_first_knearest(int zone[3])
{
    axes[2][0]=double(zone[0]); axes[2][1]=double(zone[1]); axes[2][2]=double(zone[2]);
    vector_normalize(axes[2], axes[2]);
    vector_zero(axes[2], axes[2]);
    vector_copy(axes[0], knearest_1->K);
    vector_normalize(axes[0], axes[0]);
    vector_zero(axes[0], axes[0]);
    vector_cross(axes[1], axes[0], axes[2]);
    vector_normalize(axes[1], axes[1]);
    vector_zero(axes[1], axes[1]);
    vector_cross(axes[0], axes[1], axes[2]);
    vector_normalize(axes[0], axes[0]);
    vector_zero(axes[0], axes[0]);
}

void KED::rotate(double x[3], double y[3])
{
    vector_copy(axes[0], x);
    vector_normalize(axes[0], axes[0]);
    vector_zero(axes[0], axes[0]);
    vector_copy(axes[1], y);
    vector_normalize(axes[1], axes[1]);
    vector_zero(axes[1], axes[1]);
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f]", x[0], x[1], x[2], y[0], y[1], y[2]);
        exit(1);
    }
}

void KED::ked(char *ked_path)
{
    if(mpi_rank!=0) return;
    FILE *fp=nullptr;
    fp=fopen(ked_path,"w");
    fprintf(fp, "# N_1\tN_2\tN_3\tK_1\tK_2\tK_3\tx\ty\tz\tintensity\tintensity_norm (%d points, rotated by x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], and z-[%.8f %.8f %.8f])\n", numk-1,
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double *pos_x=nullptr, *pos_y=nullptr, *intensity=nullptr;
    callocate(&pos_x, numk-1, 0.0); 
    callocate(&pos_y, numk-1, 0.0); 
    callocate(&intensity, numk-1, 0.0);
    KED_KNODE *ktemp=khead->next;
    double constn=100.0/intensity_max;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        double intensity_norm=constn*ktemp->intensity;
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), ktemp->K[0], ktemp->K[1], ktemp->K[2], xyz[0], xyz[1], xyz[2], ktemp->intensity, intensity_norm);
        fflush(fp);
        pos_x[i-1]=xyz[0]; pos_y[i-1]=xyz[1]; intensity[i-1]=intensity_norm;
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", ked_path);

    char png_path[strlen(ked_path)+5];
    strcpy(png_path, ked_path); strcat(png_path, ".png");
    img(png_path, pos_x, pos_y, intensity, numk-1, Kmagnitude_max);
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void KED::ked(char *ked_path, double sigma, double dx)
{
    if(mpi_rank!=0) return;
    int    nbin=2*round(Kmagnitude_max/dx)+1;
    int    nbin_half=nbin/2;
    double *pos_x=nullptr, **intensity=nullptr;
    callocate(&pos_x, nbin, 0.0);
    callocate_2d(&intensity, nbin, nbin, 0.0);
    for(int i=0;i<=nbin_half;i++){
        pos_x[nbin_half+i]=double(i)*dx; 
        pos_x[nbin_half-i]=-double(i)*dx;
    }

    double **intensity_c=nullptr; 
    callocate_2d(&intensity_c, nbin, nbin, 0.0);
    KED_KNODE *ktemp=khead->next;
    for(int i=1;i<numk&&ktemp!=nullptr;i++){
        double xyz[3]; vector_rotate(xyz, axes, ktemp->K);
        int m=round(xyz[0]/dx)+nbin_half;
        int n=round(xyz[1]/dx)+nbin_half;
        gaussian(intensity_c, pos_x, pos_x, nbin, ktemp->intensity, pos_x[m], pos_x[n], sigma);
        for(int j=0;j<nbin;j++){
            for(int k=0;k<nbin;k++){
                intensity[j][k]+=intensity_c[j][k];
                intensity_c[j][k]=0.0;
            }
        }
        ktemp=ktemp->next;
    }
    deallocate_2d(intensity_c, nbin);

    double imax=0.0, imin=1.0e8;
    for(int i=0;i<nbin;i++){
        for(int j=0;j<nbin;j++){
            if(imax<intensity[i][j]) imax=intensity[i][j];
            if(imin>intensity[i][j]) imin=intensity[i][j];
        }
    }
    printf("[INFO] Number of profiled diffraction intensity: %d\n", nbin*nbin);
    printf("[INFO] Range of profiled diffraction intensity: %.8f %.8f\n", imin, imax);

    FILE *fp=nullptr;
    fp=fopen(ked_path,"w");
    fprintf(fp, "# x\ty\tintensity\tintensity_norm (%d points, rotated by x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], and z-[%.8f %.8f %.8f])\n", nbin*nbin,
            axes[0][0], axes[0][1], axes[0][2], axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
    double constn=100.0/imax;
    for(int i=0;i<nbin;i++){
        for(int j=0;j<nbin;j++){
            fprintf(fp, "%.8f\t%.8f\t%.8f\t%.8f\n", pos_x[j], pos_x[i], intensity[i][j], constn*intensity[i][j]);
            fflush(fp);
        }
    }
    fclose(fp);
    printf("[INFO] Information for diffraction pattern stored in %s\n", ked_path);
    
    char png_path[strlen(ked_path)+5];
    strcpy(png_path, ked_path); strcat(png_path, ".png");
    image_array(png_path, intensity, imax, imin, nbin, nbin, 'w');
    printf("[INFO] Image for diffraction pattern stored in %s\n", png_path);
}

void KED::img(char *png_path, double *x, double *y, double *value, int num, double limit)
{
    double height=6.0, width=6.0;
    int tick_max=int(limit);
    int n_major_tick=2*tick_max+1;
    double *major_ticks; mallocate(&major_ticks, n_major_tick);
    for(int i=0;i<n_major_tick;i++){
        major_ticks[i]=double(-tick_max+i);
    }
    GRAPH graph(width, height, 300);
    graph.set_xlim(-limit, limit);
    graph.set_ylim(-limit, limit);
    graph.set_xticks(major_ticks, n_major_tick);
    graph.set_yticks(major_ticks, n_major_tick);
    graph.set_tick_in(false);
    graph.scatter(x, y, value, num);
    graph.draw(png_path);
}

KED::KED(char *ked3_path)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if(mpi_rank!=0) return;
    FILE *fp=fopen(ked3_path, "r");
    if(fp==NULL){
        printf("[ERROR] Unable to open file %s\n", ked3_path);
        exit(1);
    }
    fseek(fp, 0, SEEK_SET);

    char readbuff[MAX_LENTH_IN_NAME];
    int  keyi=0;
    for(int i=0;(keyi<2)&&i<MAX_WORD_NUMBER;i++){
        fscanf(fp, "%s", readbuff);
        if('#'==readbuff[0]){
            while(('\n'!=fgetc(fp))&&(!feof(fp)));
            continue;
        }
        if(0==strcmp(readbuff, "ELECTRON_WAVELENTH")){
            fscanf(fp, "%lf", &lambda);
            keyi++;
        }
        if(0==strcmp(readbuff, "DIFFRACTION_INTENSITY")){
            int num;
            fscanf(fp, "%d", &num);
            double hkl[3], K[3];
            double Kmagnitude, intensity;
            for(int i=0;i<num;i++){
                fscanf(fp, "%lf", &hkl[0]);
                fscanf(fp, "%lf", &hkl[1]);
                fscanf(fp, "%lf", &hkl[2]);
                fscanf(fp, "%lf", &K[0]);
                fscanf(fp, "%lf", &K[1]);
                fscanf(fp, "%lf", &K[2]);
                fscanf(fp, "%lf", &Kmagnitude);
                fscanf(fp, "%lf", &intensity);
                add_k_node(hkl, K, Kmagnitude, intensity);
                if(Kmagnitude_max<Kmagnitude) Kmagnitude_max=Kmagnitude;
            }
            keyi++;
        }
    }
    if(keyi<2){
        printf("[ERROR] Unable to read file %s\n", ked3_path);
        exit(1);
    }
}

void KED::split_k_node()
{
    struct NODE{double h,k,l; double K1, K2, K3; double Kmagnitude; double intensity;};
    MPI_Datatype MPI_NODE;
    int blocklengths[4] = {3, 3, 1, 1};
    MPI_Aint offsets[4] = {offsetof(NODE, h), offsetof(NODE, K1), offsetof(NODE, Kmagnitude), offsetof(NODE, intensity)};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(4, blocklengths, offsets, types, &MPI_NODE);
    MPI_Type_commit(&MPI_NODE);

    std::vector<NODE> vec_all;
    std::vector<int> counts(mpi_size), disps(mpi_size);
    if(mpi_rank==0){
        KED_KNODE *ktemp=khead;
        while(ktemp){
            vec_all.push_back({ktemp->hkl[0],ktemp->hkl[1],ktemp->hkl[2],ktemp->K[0],ktemp->K[1],ktemp->K[2],ktemp->Kmagnitude,ktemp->intensity});
            ktemp=ktemp->next;
        }
        int num_all=vec_all.size();
        int base=numk/mpi_size, extra=numk%mpi_size;
        for(int i=0; i<mpi_size; i++)
            counts[i]=base+(i<extra?1:0);
        disps[0]=0;
        for(int i=1; i<mpi_size; i++)
            disps[i]=disps[i-1]+counts[i-1];
        free_k_node();
    }
    MPI_Bcast(counts.data(), mpi_size, MPI_INT, 0, MPI_COMM_WORLD);

    int num_local=counts[mpi_rank];
    std::vector<NODE> vec_local(num_local);
    MPI_Scatterv(vec_all.data(), counts.data(), disps.data(), MPI_NODE, vec_local.data(), num_local, MPI_NODE,0, MPI_COMM_WORLD);
    for(auto &vec:vec_all){
        double hkl[3]={vec.h, vec.k, vec.l};
        double K[3]={vec.K1, vec.K2, vec.K3};
        add_k_node(hkl, K, vec.Kmagnitude, vec.intensity);
    }
    MPI_Type_free(&MPI_NODE);
}

void KED::ked3(char *ked3_path)
{
    if(mpi_rank!=0) return;
    FILE *fp=fopen(ked3_path, "w");
    fprintf(fp, "ELECTRON_WAVELENTH\n");
    fprintf(fp, "%.8f\n", lambda);
    fprintf(fp, "DIFFRACTION_INTENSITY\n");
    fprintf(fp, "%d\n", numk);
    KED_KNODE *ktemp=khead;
    for(int i=0;i<numk&&ktemp!=nullptr;i++){
        fprintf(fp, "%d\t%d\t%d\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n", int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->Kmagnitude, ktemp->intensity);
        fflush(fp);
        ktemp=ktemp->next;
    }
    fclose(fp);
    printf("[INFO] Three-dimensional %d diffraction intensity for calculating Kikuchi pattern stored in %s\n", numk, ked3_path);
}

// void KED::vtk(char *vtk_path)
// {
//     int nump, dimension[3];
//     for(int i=0;i<3;i++){
//         dimension[i]=kmax[i]-kmin[i]+1;
//     }
//     FILE *fp=nullptr;
//     fp=fopen(vtk_path,"w");
//     fprintf(fp, "# vtk DataFile Version 3.0\n");
//     fprintf(fp, "ELECTRON DIFFRACTION INTENSITY DISTRUBUTION\n");
//     fprintf(fp, "ASCII\n");
//     fprintf(fp, "DATASET STRUCTURED_POINTS\n");
//     fprintf(fp, "DIMENSIONS %d %d %d\n", dimension[0],  dimension[1], dimension[2]);
//     fprintf(fp, "ASPECT_RATIO %g %g %g\n", spacingK[0], spacingK[1], spacingK[2]);
//     fprintf(fp, "ORIGIN %g %g %g\n", kmin[0]*spacingK[0], kmin[1]*spacingK[1], kmin[2]*spacingK[2]);
//     fprintf(fp, "POINT_DATA %d\n", dimension[0]*dimension[1]*dimension[2]);
//     fprintf(fp, "SCALARS intensity float\n");
//     fprintf(fp, "LOOKUP_TABLE default\n");
//     double ***data; 
//     callocate_3d(&data, dimension[0], dimension[1], dimension[2], -1.0);
//     KED_KNODE *ktemp=khead;
//     for(int i=0;i<numk&&ktemp!=nullptr;i++){
//         int ix=ktemp->hkl[0]-kmin[0];
//         int iy=ktemp->hkl[1]-kmin[1];
//         int iz=ktemp->hkl[2]-kmin[2];
//         data[ix][iy][iz]=ktemp->intensity;
//         ktemp=ktemp->next;
//     }
//     for(int il=0;il<dimension[2];il++){
//         for(int ik=0;ik<dimension[1];ik++){
//             for(int ih=0;ih<dimension[0];ih++){
//                 fprintf(fp, "%g\n", data[ih][ik][il]);
//                 fflush(fp);
//             }
//         }
//     }
//     deallocate_3d(data, dimension[0], dimension[1]);
//     fclose(fp);
//     printf("[INFO] Visualized data for three-dimensional kinematic electron pattern stored in %s.\n", vtk_path);
// }

KKD::KKD(KED *ked, double xaxis[3], double yaxis[3], double zaxis[3], double thickness, double ratiox, double ratioy, int npx, int npy, char *mode)
{
    mpi_rank=ked->mpi_rank; mpi_size=ked->mpi_size;
    if(mpi_rank==0){
        printf("[INFO] Starting computation of kinematic Kikuchi diffraction...\n");
        printf("[INFO] Electron wavelength [Angstrom]: %.8f\n", ked->lambda);
        printf("[INFO] Intensity at the transmission spot: %.8f\n", ked->khead->intensity);
        printf("[INFO] Total number of diffraction intensity (including intensity at the transmission spot): %d\n", ked->numk);
        printf("[INFO] Total range of diffraction intensity: %.8f %.8f\n", ked->intensity_min, ked->intensity_max);
    }
    ked->split_k_node();

    double kn=1.0/ked->lambda;
    numpx=npx; numpy=npy;
    rotate(xaxis, yaxis, zaxis);
    compute_Kikuchi_sphere_projection(ratiox, ratioy, kn, mode);
    if(mpi_rank==0){
        printf("[INFO] Kikuchi pattern has %d pixels along x-[%.8f %.8f %.8f] and %d pixels along y-[%.8f %.8f %.8f] under zone-[%.8f %.8f %.8f]\n", 
                numpx, axes[0][0], axes[0][1], axes[0][2], numpy, axes[1][0], axes[1][1], axes[1][2], axes[2][0], axes[2][1], axes[2][2]);
        printf("[INFO] Kikuchi pattern has %.8f distance [degree] along x axis and %.8f distance [degree] along y axis\n", thetax*RAD_TO_DEG, thetay*RAD_TO_DEG);
        printf("[INFO] Kikuchi pattern has %.8f distance [Angstrom-1] along x axis and %.8f distance [Angstrom-1] along y axis\n", thetax*kn, thetay*kn);
        printf("[INFO] Starting projection of diffraction intensity on the Kikuchi pattern with %d processes...\n", mpi_size);
    }

    callocate_2d(&screenI, numpy, numpx, 0.0);
    clock_t start, finish;
    start=clock();
    quick_sort(ked->khead);
    KED_KNODE *ktemp=ked->khead;
    int num=ked->numk;
    if(mpi_rank==0){
        ktemp=ktemp->next;
        num=num-1;
    }
    int my_count=0;
    bool is_count=false;
    for(int i=0;i<num&&ktemp!=nullptr;i++){
        double upper_bound=sqrt(kn*kn+ktemp->Kmagnitude*thickness/2.0);
        double lower_bound=sqrt(kn*kn-ktemp->Kmagnitude*thickness/2.0);
        // double upper_bound=ktemp->Kmagnitude/2.0+thickness/2.0;
        // double lower_bound=ktemp->Kmagnitude/2.0-thickness/2.0;
        double hkl[3]; vector_copy(hkl, ktemp->hkl);
        double intensity=ktemp->intensity;
        for(int j=0;j<numpy;j++){
            for(int k=0;k<numpx;k++){
                double d[3];
                vector_difference(d, screenK0[j][k], ktemp->K);
                double proj=vector_length(d);
                // double proj=vector_dot(screenK0[j][k], ktemp->K)/ktemp->Kmagnitude;
                if(proj<=upper_bound&&proj>=lower_bound){
                    screenI[j][k]+=intensity;
                    is_count=true;
                }
            }
        }
        ktemp=ktemp->next;
        if(is_count){
            if((0==int(hkl[0])+int(ktemp->hkl[0]))&&(0==int(hkl[1])+int(ktemp->hkl[1]))&&(0==int(hkl[2])+int(ktemp->hkl[2]))){
                for(int j=0;j<numpy;j++){
                    for(int k=0;k<numpx;k++){
                        double d[3];
                        vector_difference(d, screenK0[j][k], ktemp->K);
                        double proj=vector_length(d);
                        // double proj=vector_dot(screenK0[j][k], ktemp->K)/ktemp->Kmagnitude;
                        if(proj<=upper_bound&&proj>=lower_bound){
                            screenI[j][k]+=ktemp->intensity;
                        }
                    }
                }
                add_k_node(ktemp->hkl, ktemp->K, ktemp->Kmagnitude, ktemp->intensity, intensity);
            }
        }else if((zaxis[1]*ktemp->K[2]-zaxis[2]*ktemp->K[1])<1.0e-6&&(zaxis[2]*ktemp->K[0]-zaxis[0]*ktemp->K[2])<1.0e-6&&(zaxis[0]*ktemp->K[1]-zaxis[1]*ktemp->K[0])<1.0e-6){
            for(int j=0;j<numpy;j++){
                for(int k=0;k<numpx;k++){
                    double d[3];
                    vector_difference(d, screenK0[j][k], ktemp->K);
                    double proj=vector_length(d);
                    // double proj=vector_dot(screenK0[j][k], ktemp->K)/ktemp->Kmagnitude;
                    if(proj<=upper_bound&&proj>=lower_bound){
                        screenI[j][k]+=ktemp->intensity;
                    }
                }
            }
            add_k_node(ktemp->hkl, ktemp->K, ktemp->Kmagnitude, ktemp->intensity, intensity);
        }
        ktemp=ktemp->next;
        is_count=false;
        my_count++;
    }
    finish=clock();
    double my_time=double(finish-start)/CLOCKS_PER_SEC;
    int total_count=0;
    MPI_Reduce(&my_count, &total_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    double total_time=0.0;
    MPI_Reduce(&my_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(mpi_rank==0){
        printf("[INFO] Ending projection of diffraction intensity on the Kikuchi pattern with %d k-points and %d processes\n", total_count, mpi_size);
        printf("[INFO] Projection time [s]: %.8f.\n", total_time);
    }
    merge_screenI();
    if(mpi_rank==0){
        for(int i=0;i<numpy;i++){
            for(int j=0;j<numpx;j++){
                if(intensity_min>screenI[i][j]) intensity_min=screenI[i][j];
                if(intensity_max<screenI[i][j]) intensity_max=screenI[i][j];
            }
        }
        printf("[INFO] Range of Kikuchi intensity on Kikuchi pattern: %.8f %.8f\n", intensity_min, intensity_max);
    }
    merge_k_node();
    if(mpi_rank==0){
        printf("[INFO] Number of Kikuchi band on Kikuchi pattern: %d\n", numk);
        KKD_KNODE *ktemp=khead;
        for(int i=0;i<numk&&ktemp!=nullptr;i++){
            printf("[INFO] Kikuchi band %d: N-[%d %d %d] K-[%.8f %.8f %.8f] Kwidth-%.8f Kintensity1-%.8f Kintensity2-%.8f\n", i+1, 
            int(ktemp->hkl[0]), int(ktemp->hkl[1]), int(ktemp->hkl[2]), ktemp->K[0], ktemp->K[1], ktemp->K[2], ktemp->Kwidth, ktemp->intensity1, ktemp->intensity2);
            ktemp=ktemp->next;
        }
        printf("[INFO] Ending computation of kinematic Kikuchi diffraction\n");
    }
}

KKD::~KKD()
{
    free_screenI();
    free_k_node();
}

void KKD::rotate(double x[3], double y[3], double z[3])
{
    vector_copy(axes[0], x);
    vector_normalize(axes[0], axes[0]);
    vector_copy(axes[1], y);
    vector_normalize(axes[1], axes[1]);
    vector_copy(axes[2], z);
    vector_normalize(axes[2], axes[2]);
    for(int i=0;i<3;i++){
        vector_zero(axes[i], axes[i]);
    }
    if(fabs(vector_dot(axes[0], axes[1]))>1.0e-6||fabs(vector_dot(axes[1], axes[2]))>1.0e-6||fabs(vector_dot(axes[0], axes[2]))>1.0e-6){
        printf("[ERROR] The orthogonality condition is not satisfied with x-[%.8f %.8f %.8f], y-[%.8f %.8f %.8f], z-[%.8f %.8f %.8f]", 
                x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2]);
        exit(1);
    }
}

void KKD::compute_Kikuchi_sphere_projection(double ratiox, double ratioy, double kn, char *mode)
{
    callocate_3d(&screenK0, numpy, numpx, 3, 0.0);
    int impx=numpx/2, impy=numpy/2;
    int err;
    if(0==strcmp(mode, "stereo")){
        for(int i=0;i<numpy;i++){
            for(int j=0;j<numpx;j++){
                double xy[2]={double(i-impy)/double(impy)*ratioy, double(j-impx)/double(impx)*ratiox};
                double xyz[3];
                compute_sphere_from_stereographic_projection(xyz, err, xy);
                if(0==err){
                    vector_transform(xyz, xyz, axes);
                    vector_copy(screenK0[i][j], xyz);
                }
            }
        }
    }else if(0==strcmp(mode, "ortho")){
        for(int i=0;i<numpy;i++){
            for(int j=0;j<numpx;j++){
                double xy[2]={double(i-impy)/double(impy)*ratioy, double(j-impx)/double(impx)*ratiox};
                double xyz[3];
                compute_sphere_from_orthographic_projection(xyz, err, xy);
                if(0==err){
                    vector_transform(xyz, xyz, axes);
                    vector_copy(screenK0[i][j], xyz);
                }
            }
        }
    }else{
        printf("[ERROR] Unrecognized projection mode %s\n", mode);
        exit(1);
    }
    thetax=acos(vector_dot(screenK0[impy][0], screenK0[impy][numpx-1]));
    thetay=acos(vector_dot(screenK0[0][impx], screenK0[numpy-1][impx]));
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            vector_constant(screenK0[i][j], kn, screenK0[i][j]);
        }
    }
}

void KKD::merge_screenI()
{
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpi_rank==0){
        double *row_buffer; mallocate(&row_buffer, numpx);
        for(int row=0; row<numpy; row++){
            if(row==0){
                MPI_Reduce(MPI_IN_PLACE, screenI[row], numpx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            } else {
                MPI_Reduce(screenI[row], row_buffer, numpx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                for(int col=0; col<numpx; col++){
                    screenI[row][col]=row_buffer[col];
                }
            }
        }
        deallocate(row_buffer); 
    }else{
        for(int row=0; row<numpy; row++){
            MPI_Reduce(screenI[row], NULL, numpx, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        free_screenI();
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void KKD::free_screenI()
{
    if(0!=numpx&&0!=numpy){
        deallocate_2d(screenI, numpy);
        deallocate_3d(screenK0, numpy, numpx);
    }
    numpy=numpx=0;
}

void KKD::add_k_node(double hkl[3], double K[3], double Kwidth, double intensity1, double intensity2)
{
    if(khead==nullptr&&ktail==nullptr){
        khead=ktail=new KKD_KNODE;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kwidth=Kwidth;
        ktail->intensity1=intensity1;
        ktail->intensity2=intensity2;
        numk++;
    }else{
        ktail->next=new KKD_KNODE;
        ktail=ktail->next;
        vector_copy(ktail->hkl, hkl);
        vector_copy(ktail->K, K);
        vector_zero(ktail->K, ktail->K);
        ktail->Kwidth=Kwidth;
        ktail->intensity1=intensity1;
        ktail->intensity2=intensity2;
        numk++;
    }
}

void KKD::merge_k_node()
{
    struct NODE{double h,k,l; double K1, K2, K3; double Kwidth; double intensity1, intensity2;};
    MPI_Datatype MPI_NODE;
    int blocklengths[4] = {3, 3, 1, 2};
    MPI_Aint offsets[4] = {offsetof(NODE, h), offsetof(NODE, K1), offsetof(NODE, Kwidth), offsetof(NODE, intensity1)};
    MPI_Datatype types[4] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Type_create_struct(4, blocklengths, offsets, types, &MPI_NODE);
    MPI_Type_commit(&MPI_NODE);

    int numk_all=0;
    MPI_Reduce(&numk, &numk_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    std::vector<int> numks(mpi_size);
    MPI_Gather(&numk,1,MPI_INT,numks.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    std::vector<NODE> vec;
    KKD_KNODE *ktemp=khead;
    while(ktemp){
        vec.push_back({ktemp->hkl[0],ktemp->hkl[1],ktemp->hkl[2],ktemp->K[0],ktemp->K[1],ktemp->K[2],ktemp->Kwidth,ktemp->intensity1,ktemp->intensity2});
        ktemp=ktemp->next;
    }

    std::vector<int> disps(mpi_size);
    std::vector<NODE> vec_all;
    if(mpi_rank==0){
        vec_all.resize(numk_all);
        disps[0]=0;
        for(int i=1; i<mpi_size; ++i){
            disps[i]=disps[i-1]+numks[i-1];
        }
    }
    MPI_Gatherv(vec.data(), numk, MPI_NODE, vec_all.data(), numks.data(), disps.data(), MPI_NODE, 0, MPI_COMM_WORLD);
    free_k_node();
    
    if(mpi_rank==0){
        for(auto &vec : vec_all){
            double hkl[3]={vec.h,vec.k,vec.l};
            double K[3]={vec.K1,vec.K2,vec.K3};
            add_k_node(hkl, K, vec.Kwidth, vec.intensity1, vec.intensity2);
        }
    }
    MPI_Type_free(&MPI_NODE);
}

void KKD::free_k_node()
{
    KKD_KNODE *cur=khead;
    while(cur!=nullptr){
        KKD_KNODE *temp=cur;
        cur=cur->next;
        delete temp;
    }
    khead=ktail=nullptr;
    numk=0;  
}

void KKD::kkd(char* kkd_path, char background)
{
    if(mpi_rank!=0) return;
    FILE *fp=nullptr;
    fp=fopen(kkd_path,"w");
    fprintf(fp, "KIKUCHI_IMAGE_SIZE\n");
    fprintf(fp, "%d %d\n", numpx, numpy);
    fprintf(fp, "KIKUCHI_IMAGE_VALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenI[i][j]);
            fflush(fp);
        }
    }
    printf("[INFO] Information for Kikuchi pattern stored in %s\n", kkd_path);
    char png_path[strlen(kkd_path)+5];
    strcpy(png_path, kkd_path); strcat(png_path, ".png");
    image_array(png_path, screenI, intensity_max, intensity_min, numpy, numpx, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
}

void KKD::kkd(char* kkd_path, double vmax, double vmin, char background)
{
    if(mpi_rank!=0) return;
    FILE *fp=nullptr;
    fp=fopen(kkd_path,"w");
    fprintf(fp, "KIKUCHI_IMAGE_SIZE\n");
    fprintf(fp, "%d %d\n", numpx, numpy);
    fprintf(fp, "KIKUCHI_IMAGE_VALUE\n");
    fflush(fp);
    for(int i=0;i<numpy;i++){
        for(int j=0;j<numpx;j++){
            fprintf(fp, "%.8f\n", screenI[i][j]);
            fflush(fp);
        }
    }
    printf("[INFO] Information for Kikuchi pattern stored in %s\n", kkd_path);
    char png_path[strlen(kkd_path)+5];
    strcpy(png_path, kkd_path); strcat(png_path, ".png");
    image_array(png_path, screenI, vmax, vmin, numpy, numpx, background);
    printf("[INFO] Image for Kikuchi pattern stored in %s\n", png_path);
}