#ifndef __AAVDP_XRD_H__
#define __AAVDP_XRD_H__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <mpi.h>
#include "../MATH/MATH.h"
#include "../MATH/GRAPH.h"
#include "../MODEL/MODEL.h"

#define XRD_INTENSITY_LIMIT 1.0e-6
struct XRD_KNODE{
    int    hkl[3];
    double theta;//radian
    double intensity;
    double deviation;
    int    multiplicity;
    XRD_KNODE *next=nullptr;
};

extern void copy_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2);
extern void swap_knode_data(XRD_KNODE *knode1, XRD_KNODE *knode2);
extern void quick_sort(XRD_KNODE *kstart, XRD_KNODE *kend);

#define SCHERRER_CONST 0.90


//x-ray diffraction
class XRD
{
public:
    int    mpi_rank=0, mpi_size=1;
    int    numk=0;
    XRD_KNODE  *khead=nullptr;
    XRD_KNODE  *ktail=nullptr;
    double minTheta, maxTheta;
    double intensity_min=1.0e8, intensity_max=0.0;
    XRD(MODEL *model, double min2Theta, double max2Theta, double bin2Theta, double threshold, double spacing[3], bool is_spacing_auto);
    ~XRD();
    void   xrd(char *xrd_path);
    void   xrd(char *xrd_path, double NA, double NB, double U, double V, double W, double bin2Theta);
    // void   xrd(char *xrd_path, double mixing_param, double scherrer_lambda, double scherrer_diameter, double bin2Theta);
    // void   xrd(char *xrd_path, double mixing_param, double FWHM, double bin2Theta);
private:
    void   add_k_node(int hkl[3], double theta, double intensity, double deviation, int multiplicity);
    void   merge_k_node();
    void   free_k_node();
    void   filter_diffraction_intensity(double threshold);
    // void   unique_diffraction_intensity();
    void   histogram_diffraction_intensity(double bin2Theta);
    void   img(char *png_path, double *x, double *y, int num, double xmin, double xmax, char mode);
};

#endif