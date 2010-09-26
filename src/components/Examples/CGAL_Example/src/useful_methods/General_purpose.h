///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2009
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef General_purpose_H
#define General_purpose_H

#include "../../../../../mepp/Polyhedron/polyhedron.h"
#include "../triangle_tools/Triangle_shape_analysis.h"

#include "time.h"
#include <string>
#include <vector>
///////////////////////////////////////////////////////////////////////////////////////////////////////
// for files...
static inline void open_file(const char* file_name, FILE *&pFile, char* open_option)
{
    pFile = fopen(file_name,open_option);
    if (!pFile)
    {
        cerr << "Cannot open the file: " << file_name << endl;
        exit(-1);
    }
}

static void load_testing_results_TT_TF_FF_FT(   const char* file_name,
                                                unsigned int* nb_TT_G,
                                                unsigned int* nb_TF_G,
                                                unsigned int* nb_FF_G,
                                                unsigned int* nb_FT_G,
                                                unsigned int nb_lines
                                                )
{
    FILE *pFile = fopen(file_name,"r");

    unsigned int nb_line=0;
    if (pFile)
    {

        while(fscanf(pFile, "%d %d %d %d", &nb_TT_G[nb_line],
                                           &nb_TF_G[nb_line],
                                           &nb_FF_G[nb_line],
                                           &nb_FT_G[nb_line]) != EOF)
                                           nb_line++;

        fclose(pFile);
    }
    else
    {
        while(nb_line<nb_lines)
        {
            nb_TT_G[nb_line] = nb_TF_G[nb_line] = nb_FF_G[nb_line] = nb_FT_G[nb_line] = 0;
            nb_line++;
        }
    }

    cout << "nb_lines in " << file_name << ": " << nb_line << endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
// displaying elapsed time
static inline void print_elapsed_time(const long& start_time_s)
{
    unsigned long time_s = clock()/CLOCKS_PER_SEC - start_time_s;
	double time_m = time_s/60;
	time_s = time_s - time_m*60;
	cout << " Time = " << time_m << " min " << time_s << " secondes" << endl;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
// vectors
static inline void print_vector(const Vector& V)
{
	cout << " (" << V.x() << ", " << V.y() << ", " << V.z() << ")" << endl;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////
// image methods (DEM...)

static inline void random_8bits_image(vector <unsigned int>& DEM, unsigned int& nbline, unsigned int& nbcolumn)
{
    DEM.clear();
    srand((unsigned)time(NULL));
    nbline = nbcolumn = 512;
    for(unsigned int i=0; i<nbline; ++i)
        for(unsigned int j=0; j<nbcolumn; ++j)
        {
            DEM.push_back((unsigned int)(rand()/double(RAND_MAX) * 255));
        }
}

static inline void load_DEM_textfile(const char* text_file_name, vector <unsigned int>& DEM, unsigned int& nbline, unsigned int& nbcolumn)
{
    DEM.clear();

    ifstream flux(text_file_name, ios::in);
    // open the file:
    if(!flux)
    {
        cout << "read_triangle_indices: Cannot open the file: "<< text_file_name << "!!!" << endl;
        //exit(-1); // we do not exit here
        return;
    }
    string line;
    unsigned int cpt=0;
    while ( getline(flux, line) )
    { // for each triangle
        istringstream iss(line);
        double i;
        iss >> i;
        //cout << i << endl;
        if(cpt>1)
        {
            if(cpt>=nbline*nbcolumn+2) break; // we do not take into account the remaining values if there are
            DEM.push_back(i);
        }
        else if(cpt==0)
            nbline = (unsigned int)i;
        else nbcolumn = (unsigned int)i;

        cpt++;
    }
    //cout << "nbline = " << nbline << " nbcolumn = " << nbcolumn << " nbline*nbcolumn = " << nbline*nbcolumn << " DEM.size() = " << DEM.size() << endl;
    assert(DEM.size()==nbline*nbcolumn);
    flux.close();
}

/*
Currently the following file formats are supported:

    * Windows bitmaps - BMP, DIB
    * JPEG files - JPEG, JPG, JPE
    * Portable Network Graphics - PNG
    * Portable image format - PBM, PGM, PPM
    * Sun rasters - SR, RAS
    * TIFF files - TIFF, TIF

static inline void load_DEM(   const char* image_file_name,
                        vector <unsigned int>& DEM,
                        unsigned int& nbline,
                        unsigned int& nbcolumn)
{
    //IplImage *buf = cvCreateImage(cvSize(imWidth, imHeight), IPL_DEPTH_8U, 1);

    // the following line needs highgui
    IplImage * image;// = cvLoadImage(image_file_name, 0); // we load the image in grey scale level!! 1 chanel image

    if(!image)
    {
        cout << "load_DEM: cannot find image file." << endl;
        exit(-1);
    }

    nbline = image->height;
    nbcolumn = image->width;

    unsigned char* ptdata = (unsigned char*)image->imageData;
    unsigned int cpt=0;
    while(cpt<nbline*nbcolumn)
    {
        DEM.push_back(*ptdata);
        ptdata++;
        cpt++;
    }
    cvReleaseImage(&image);
}
*/
static inline double linear_interpolation_image(   unsigned int x1, unsigned int y1, // pixel left (x1,y1)
                                            unsigned int x2, unsigned int y2, // pixel right (x2,y2)
                                            double val_left,
                                            double val_right,
                                            double x, double y        // continuous position to interpolate
                                            )
{
    assert((x1 < x2 && y1==y2) || (y1 < y2 && x1==x2));
    if(y1==y2)
    { // we interpolate linearly along Ox
        return (x2-x)/(x2-x1)*val_left + (x-x1)/(x2-x1)*val_right;
    }
    else //if(x1==x2)
    { // we interpolate linearly along Oy
        return (y2-y)/(y2-y1)*val_left + (y-y1)/(y2-y1)*val_right;
    }
};

static inline double bilinear_interpolation_image( unsigned int x1, unsigned int y1, // pixel bottom left corner (x1,y1)
                                            unsigned int x2, unsigned int y2, // pixel top right corner (x2,y2)
                                            double val_bottom_left,
                                            double val_bottom_right,
                                            double val_top_right,
                                            double val_top_left,
                                            double x, double y        // continuous position to interpolate
                                            )
{
    assert(x1 < x2);
    assert(y1 < y2);
    // we first interpolate linearly along Ox
    double f12 = (x2-x)/(x2-x1)*val_bottom_left + (x-x1)/(x2-x1)*val_bottom_right;
    double f34 = (x2-x)/(x2-x1)*val_top_left + (x-x1)/(x2-x1)*val_top_right;
    //cout << "bilinear_interpolation_image: val_bottom_left = " << val_bottom_left << " val_bottom_right = " << val_bottom_right << endl;
    //cout << "bilinear_interpolation_image: val_top_left = " << val_top_left << " val_top_right = " << val_top_right << endl;
    //cout << "bilinear_interpolation_image: f12 = " << f12 << " f34 = " << f34 << endl;
    // we then interpolate along Oy
    return (y2-y)/(y2-y1)*f12 + (y-y1)/(y2-y1)*f34;
};


static inline double get_image_value(  const vector <unsigned int>& DEM, unsigned int nbline, unsigned int nbcolumn,
                                unsigned int x, unsigned int y) // 0x <==> column; 0y <==> lines
{
    if(DEM.empty())
    {
        cout << "get_image_value: ERROR: the image is empty." << endl;
        exit(-1);
    }
    assert(nbline*nbcolumn == DEM.size());
    //assert(x>=0); assert(x<nbcolumn);
    clamp(x, 0, nbcolumn);
    //assert(y>=0); assert(y<nbline);
    clamp(y, 0, nbline);

    return DEM[y*nbcolumn+x];
};

static inline double get_image_value(  const vector <unsigned int>& DEM, unsigned int nbline, unsigned int nbcolumn,
                                double x, double y) // 0x <==> column; 0y <==> lines
{ // x number_type_delaunayresents column index, y represents line index (we think of the origin at (0,0) as in the cartesian system)
    if(DEM.empty())
    {
        cout << "get_image_value: ERROR: the image is empty." << endl;
        exit(-1);
    }
    assert(nbline*nbcolumn == DEM.size());

    // we infinitely extend the image out of its boundaries by the line boundary value:
    clamp(x, 0.0, double(nbcolumn-1));
    clamp(y, 0.0, double(nbline-1));
    //assert(x>=0.0); assert(x<nbcolumn);
    //assert(y>=0.0); assert(y<nbline);

    unsigned int x1=floor(x), y1=floor(y), x2=ceil(x), y2=ceil(y);

    if(x1==x2 && y1==y2)
    { // case where x and y are integer (no interpolation to do)
        return get_image_value(DEM, nbline, nbcolumn, x1, y1);
    }
    else if(x1==x2 || y1==y2)
    {
         return linear_interpolation_image( x1, y1, x2, y2,
                                            get_image_value(DEM, nbline, nbcolumn, x1, y1),
                                            get_image_value(DEM, nbline, nbcolumn, x2, y2), x, y);
    }

    //cout << "x = " << x << " y = " << y << endl;
    //cout << "x1 = " << x1 << " x2 = " << x2 << " y1 = " << y1 << " y2 = " << y2 << endl;
    double  val_bottom_left = get_image_value(DEM, nbline, nbcolumn, x1, y1),
            val_bottom_right = get_image_value(DEM, nbline, nbcolumn, x2, y1),
            val_top_right = get_image_value(DEM, nbline, nbcolumn, x2, y2),
            val_top_left = get_image_value(DEM, nbline, nbcolumn, x1, y2);

    return bilinear_interpolation_image(x1, y1, x2, y2, val_bottom_left, val_bottom_right, val_top_right, val_top_left, x, y);
};


static inline double get_interpolation_triangle_abs_error( const Point3d& a, const Point3d& b, const Point3d& c, // triangle in cartesian frame
                                                        const vector <unsigned int>& DEM, unsigned int nbline, unsigned int nbcolumn, // image
                                                        const Point3d& Sample) // (Sample.x(), Sample.y()) is a point in the image frame
{
    Point3d SampleCartesianFrame(Sample.x(), double(nbline)-1.0-Sample.y(),0);
    double biltri = bilinear_interpolation_triangle(a, b, c, SampleCartesianFrame,
                                                    get_image_value(DEM, nbline, nbcolumn, a.x(), double(nbline)-1.0-a.y()),
                                                    get_image_value(DEM, nbline, nbcolumn, b.x(), double(nbline)-1.0-b.y()),
                                                    get_image_value(DEM, nbline, nbcolumn, c.x(), double(nbline)-1.0-c.y()));
    //cout << "biltri = " << biltri << endl;
    double valimage = get_image_value(DEM, nbline, nbcolumn, Sample.x(), Sample.y());
    //cout << "valimage = " << valimage << endl;
    return abs(valimage-biltri);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    DEM & integral over triangle
*/
static inline double integral_over_triangle(const Point3d& a, const Point3d& b, const Point3d& c, // in cartesian frame
                                            const vector <unsigned int>& DEM, unsigned int nbline, unsigned int nbcolumn)
{ // the 3 points must be distinct and non collinear and [ab] must be the triangle base
    double sum = 0.0; // integral init
    unsigned int nbelements = 0;
    double h = 0.01; // integral precision for a square (we want to divide the triangle into small squares = pixels)
                      // with 0.001 precision I have noticed a maximum of 2% of difference between this method and the exact area
                      // and most of the time the precision is much better (0.5%)
    double h2=h*h;
    Vector Vab(b-a), Vac(c-a), Vbc(c-b); // non-normalized vectors
    //cout << "a = (" << a.x() << "," << a.y() << ") b = (" << b.x() << "," << b.y() << ") c = (" << c.x() << "," << c.y() << ")" << endl;
    double lab = std::sqrt(Vab*Vab), lac = std::sqrt(Vac*Vac), lbc = std::sqrt(Vbc*Vbc);
    Vector Vabn = Vab/lab, Vabn_ortho; // normalized vectors
    Vabn_ortho = Vector(-Vabn.y(),Vabn.x(),0); // because we are in 2D!!
    double interm = Vac*Vabn;
    //cout << "interm = " << interm << " lab = " << lab << endl;
    double dab;
    for(dab=h; dab<interm; dab+=h)
    { // along the side ab... a+dab.Vab
        // we now go in the perpendicular direction
        double dac, H;
        if(Vab*Vac !=0.0)
        {
            dac = (dab * lab*lac)/(Vab*Vac);
            H = std::sqrt(dac*dac-dab*dab);
            //cout << "(Vac*Vab) = " << (Vac*Vab) << " dac = " << dac << " dac*dac-dab*dab = " << dac*dac-dab*dab << " H = " << H << endl;
        }
        else
        {
            cout << "integral_over_triangle: some check should be done." << endl;
            H = lac;
        }

        for(double high=h; high<=H; high+=h)
        { // a + dab.Vab + high.Vab_ortho
            Point3d currentSampleCartesianFrame(a+dab*Vabn+(high-h)*Vabn_ortho);
            Point3d currentSampleImageFrame(currentSampleCartesianFrame.x(),
                                            double(nbline)-1.0-currentSampleCartesianFrame.y(),
                                            0);
            if(isInTriangle(a, b, c, currentSampleCartesianFrame))
            {
                // volume computation
                sum += h2 * get_interpolation_triangle_abs_error(a, b, c, DEM, nbline, nbcolumn, currentSampleImageFrame);
                nbelements++;
            }

            //cout << "sum = " << sum << endl;
        }
    }

    for(; dab<=lab; dab+=h)
    { // along the side ab...
        // we now go in the perpendicular direction
        double dbc, H;
        if(Vbc*-Vab!=0)
        {
            dbc = (lab-dab)*lbc*lab/(Vbc*-Vab);
            H=std::sqrt(dbc*dbc-(lab-dab)*(lab-dab));
            //cout << " dbc = " << dbc << " dbc*dbc-(lab-dab)*(lab-dab) = " << dbc*dbc-(lab-dab)*(lab-dab) << " H = " << H << endl;
        }
        else
        {
            cout << "integral_over_triangle: some check should be done." << endl;
            H = lbc;
        }

        for(double high=h; high<=H; high+=h)
        {// a + dab.Vab + high.Vab_ortho
            Point3d currentSampleCartesianFrame(a+dab*Vabn+(high-h)*Vabn_ortho);
            Point3d currentSampleImageFrame(currentSampleCartesianFrame.x(),
                                            double(nbline)-1.0-currentSampleCartesianFrame.y(),
                                            0);
            if(isInTriangle(a, b, c, currentSampleCartesianFrame))
            {
                sum += h2 * get_interpolation_triangle_abs_error(a, b, c, DEM, nbline, nbcolumn, currentSampleImageFrame);
                nbelements++;
            }

            //cout << "sum = " << sum << endl;
        }

    }

    //cout << "sum = " << sum << " triangleArea(a, b, c) = " << triangleArea(a, b, c) << endl;
    if(sum<1e-8) sum = 1e8;
    //sum = sum/(nbelements); // to avoid to give to much importance to one big error
    return sum/triangleArea(a, b, c); // we normalize the result by the triangle area to not be sensitive to the scale
};


static inline double integral_over_triangle_choosing_basis(const Point3d& a, const Point3d& b, const Point3d& c,
                                                    const vector <unsigned int>& DEM, unsigned int nbline, unsigned int nbcolumn)
{
    Vector Vab(b-a), Vac(c-a), Vbc(c-b);
    double lab2 = Vab*Vab, lac2 = Vac*Vac, lbc2 = Vbc*Vbc;
    if(lab2>=lac2 && lab2>=lbc2)
    {
        return integral_over_triangle(a, b, c, DEM, nbline, nbcolumn);
    }
    else if(lac2>=lab2 && lac2>=lbc2)
    {
        return integral_over_triangle(a, c, b, DEM, nbline, nbcolumn);
    }
    else
    {
        return integral_over_triangle(b, c, a, DEM, nbline, nbcolumn);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
// clip to sphere
static inline void clipToSphere(double sphereRadius, Vector& VtoClip)
{
    double n2 = VtoClip*VtoClip;

    if(n2!=0.0)
    {
        VtoClip = VtoClip*0.999999*sphereRadius/std::sqrt(n2); // 0.999999 to be inside the sphere
    }
}

static inline void clipToSphere(const Point3d& sphereCenter, double sphereRadius, Vector& VtoClip, Point3d& PonSphere)
{
    clipToSphere(sphereRadius, VtoClip);
    PonSphere = sphereCenter+VtoClip;
}

static inline void clipToSphere(const Point3d& sphereCenter, double sphereRadius, Vertex_handle & pVertex)
{
    Vector VtoClip = pVertex->point() - sphereCenter;
    clipToSphere(sphereCenter, sphereRadius, VtoClip, pVertex->point());
}

static inline bool isInSphere(const Point3d& sphereCenter, double sphereRadius, const Vertex_handle & pVertex)
{
    return (pVertex->point()-sphereCenter)*(pVertex->point()-sphereCenter) < sphereRadius*sphereRadius;

}

static inline bool isInSphere(const Point3d& sphereCenter, double sphereRadius, const Halfedge_handle & h)
{
    return  isInSphere(sphereCenter, sphereRadius, h->vertex()) &&
            isInSphere(sphereCenter, sphereRadius, h->prev()->vertex());

}

static inline bool isInSphere(const Point3d& sphereCenter, double sphereRadius, const Facet_handle & pFacet)
{
    Halfedge_around_facet_circulator h = pFacet->facet_begin();
    Halfedge_around_facet_circulator he = h;

    CGAL_For_all(h,he)
    {
        if( !isInSphere(sphereCenter, sphereRadius, h) ) return false;
    }

    return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
//
static bool isInside(const Vertex_handle pVertex, const vector< Vertex_handle >& points)
{
    if(pVertex==Vertex_handle()) return false;

    vector< Vertex_handle >::const_iterator it(points.begin()), ite(points.end());
    for(;it!=ite;++it)
    {
        if(*it==pVertex) return true;
    }
    return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// center-surround mechanism

// topological K-neighborhood
// can also take into account a radius to avoid taking far away points
static void N_topo( const Vertex_handle pVertex,
                    unsigned int K,
                    vector< Vertex_handle >& neighbors,
                    double sphereRadius = 1e8)
{
    neighbors.clear();

    if(pVertex==Vertex_handle()) return;

    std::set< void* > treated_Vertex;

    neighbors.push_back(pVertex);

    treated_Vertex.insert(&(*pVertex));

    if(K==0) return;

    Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
    Halfedge_around_vertex_circulator he = h;
    CGAL_For_all(h,he)
    {
        if(isInSphere(pVertex->point(), sphereRadius, h->prev()->vertex()))
        {
            neighbors.push_back(h->prev()->vertex());
            treated_Vertex.insert(&(*h->prev()->vertex()));
        }
    }

    if(K==1) return;

    h = pVertex->vertex_begin();
    he = h;
    CGAL_For_all(h,he)
    {
        Halfedge_around_vertex_circulator h2 = h->prev()->vertex()->vertex_begin();
        Halfedge_around_vertex_circulator he2 = h2;
        CGAL_For_all(h2,he2)
        {
            if( isInSphere(pVertex->point(), sphereRadius, h2->prev()->vertex()) &&
                treated_Vertex.find(&(*h2->prev()->vertex()))==treated_Vertex.end()) // if not found we add it!
            {
                neighbors.push_back(h2->prev()->vertex());
                treated_Vertex.insert(&(*h2->prev()->vertex()));
            }
        }
    }

    if(K==2) return;

    h = pVertex->vertex_begin();
    he = h;
    CGAL_For_all(h,he)
    {
        Halfedge_around_vertex_circulator h2 = h->prev()->vertex()->vertex_begin();
        Halfedge_around_vertex_circulator he2 = h2;
        CGAL_For_all(h2,he2)
        {
            Halfedge_around_vertex_circulator h3 = h2->prev()->vertex()->vertex_begin();
            Halfedge_around_vertex_circulator he3 = h3;
            CGAL_For_all(h3,he3)
            {
                if( isInSphere(pVertex->point(), sphereRadius, h3->prev()->vertex()) &&
                    treated_Vertex.find(&(*h3->prev()->vertex()))==treated_Vertex.end()) // if not found we add it!
                {
                    neighbors.push_back(h3->prev()->vertex());
                    treated_Vertex.insert(&(*h3->prev()->vertex()));
                }
            }
        }
    }

    if(K==3) return;

    h = pVertex->vertex_begin();
    he = h;
    CGAL_For_all(h,he)
    {
        Halfedge_around_vertex_circulator h2 = h->prev()->vertex()->vertex_begin();
        Halfedge_around_vertex_circulator he2 = h2;
        CGAL_For_all(h2,he2)
        {
            Halfedge_around_vertex_circulator h3 = h2->prev()->vertex()->vertex_begin();
            Halfedge_around_vertex_circulator he3 = h3;
            CGAL_For_all(h3,he3)
            {
                Halfedge_around_vertex_circulator h4 = h3->prev()->vertex()->vertex_begin();
                Halfedge_around_vertex_circulator he4 = h4;
                CGAL_For_all(h4,he4)
                {
                    if( isInSphere(pVertex->point(), sphereRadius, h4->prev()->vertex()) &&
                        treated_Vertex.find(&(*h4->prev()->vertex()))==treated_Vertex.end()) // if not found we add it!
                    {
                        neighbors.push_back(h4->prev()->vertex());
                        treated_Vertex.insert(&(*h4->prev()->vertex()));
                    }
                }
            }
        }
    }

    if(K>=5) cout << "Topological distance greater than 5 is not taken into account" << endl; // we do not study it for K greater than 4
}

// Gaussian-weighted average
static double G(const Vertex_handle pVertex,
                unsigned int K,
                double sigma)
{
    double sum=0.0, den=0.0;

    vector< Vertex_handle > neighbors;

    N_topo(pVertex, K, neighbors, 2.0 * sigma);

    vector< Vertex_handle >::const_iterator it(neighbors.begin()), ite(neighbors.end());
    for(;it!=ite;++it)
    {
        double tmp = exp(- ((*it)->point()-pVertex->point())*((*it)->point()-pVertex->point()) / (2.0*sigma*sigma) );

        den += tmp;

        sum += tmp * std::max<double>(abs((*it)->Kmax), abs((*it)->Kmin));
    }

    if(den) return sum/den;

    return 0.0;
}

static double similarity(const Vertex_handle pVertex, unsigned int K, double sigma)
{
    return abs(G(pVertex, K+1, sigma) - G(pVertex, K, sigma));
}

static double multipleScaleSimilarity(const Vertex_handle pVertex, double sigma)
{
    double sum=0.0;
    unsigned int i=0;
    while(i<4)
    {
        sum += similarity(pVertex, i, sigma);
        ++i;
    }

    return sum/i;
}
#endif // General_purpose_H
