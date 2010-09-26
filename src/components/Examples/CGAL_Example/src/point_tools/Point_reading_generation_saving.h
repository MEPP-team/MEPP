#ifndef Point_reading_generating_saving_H
#define Point_reading_generating_saving_H

#include "../useful_methods/Math_methods.h"

/*
    POINT SET GENERATION
*/
static void generate_points_uniform_grid(   std::vector< Point3d >& cloud_point,
                                            int nbline=20, int nbcol=20,
                                            double scale = 1,
                                            double percentage_to_keep = 100.0)
{
    //int cpt=0;
    srand((unsigned)time(NULL));
    // we go from 0 to nbline, nbcol to put a point at each pixel corner (e.g. for getting n col, we need n+1 points)
    for(int line=0; line<=nbline; ++line) // for each line (along Oy)
    {
        for(int col=0; col<=nbcol; ++col) // for each column (along Ox)
        {
            double alea = rand()/double(RAND_MAX);

            if(100.0*alea < percentage_to_keep) // we approximately keep the chosen percentage
            {
                //cpt++;
                cloud_point.push_back(Point3d(scale*double(col),scale*double(line),0));
            }
        }
    }
    //cout << "Percentage of kept candidates: " << 100*double(cpt)/((nbline+1)*(nbcol+1)) << "%" << endl;
}

static void generate_points_gaussian_distribution(std::vector< Point3d >& cloud_point, int nb = 20, double scale = 10.0)
{
    srand((unsigned)time(NULL));
    double N1, N2;
    for(int i=0; i<nb; ++i)
    {
        for(int j=0; j<nb; ++j)
        {
            gaussianNoise(N1, N2, false);
            cloud_point.push_back(Point3d(scale*double(i)*N1,scale*double(j)*N2,0));
        }
    }
}

/*
    POINT SET MODIFICATION BY NOISE ADDITION
*/
static void add_noise_to_points(std::vector< Point3d >& cloud_point,
                                double xmin, double xmax, double ymin, double ymax)
{
    std::vector< Point3d > cloud_point_inter;

    cloud_point_inter.reserve(cloud_point.size());

    unsigned int nbp=cloud_point.size();

    srand((unsigned)time(NULL));
    for(unsigned int p=0; p<nbp; ++p) // for each point
    {
        double noise = rand()/double(RAND_MAX);
        if(noise < 0.5)
            noise = 0.1*noise;
        else
            noise = 0.1*(noise-1.0);

        double x=cloud_point[p].x()+noise, y=cloud_point[p].y()+noise, z=0;

        if(x<xmin) x=xmin;
        else if(x>xmax) x=xmax;

        if(y<ymin) y=ymin;
        else if(y>ymax) y=ymax;

        cloud_point_inter.push_back(Point3d(x,y,z));
    }

    cloud_point.swap(cloud_point_inter);
}

static void add_noise_to_uniform_grid_points(   std::vector< Point3d >& cloud_point,
                                                int nbline=20, int nbcol=20 )
{
    std::vector< Point3d > cloud_point_inter;

    cloud_point_inter.reserve(cloud_point.size());

    unsigned int nb_point_with_noise=0;

    srand((unsigned)time(NULL));

    for(int line=0; line<=nbline; ++line) // for each line (along Oy) except boundaries
    {
        for(int col=0; col<=nbcol; ++col) // for each column (along Ox) except boundaries
        {
            double noise;
            if( (!line || (line==nbline)) && (!col || (col==nbcol)) ) // corner points
            {
                noise=0;
            }
            else if( !line || (line==nbline) || !col || (col==nbcol)) // since we have tested corner points=>boundary points
            {
                noise = rand()/double(RAND_MAX);
                if((!col && (line<nbline)) || (!line && (col<nbcol)) )
                {
                    noise = 0;//0.5*noise;
                }
                else if( ((col>0)&&(line==nbline)) || ((line>0) && (col==nbcol)) )
                {
                    noise = 0;//0.5*(noise-1.0);
                }
                else
                {
                    assert(false);
                }
            }
            else
            {
                noise = rand()/double(RAND_MAX);
                if(noise < 0.5)
                    noise = 0.5*noise;
                else
                    noise = 0.5*(noise-1.0);

                nb_point_with_noise++;
            }
            double x=cloud_point[line*(nbcol+1)+col].x()+noise, y=cloud_point[line*(nbcol+1)+col].y()+noise, z=0;
            clamp(x, cloud_point[0].x(), cloud_point[nbcol].x());
            clamp(y, cloud_point[0].y(), cloud_point[(nbline)*(nbcol+1)].y());
            cloud_point_inter.push_back(Point3d(x,y,z));
        }
    }

    cout << "add_noise_to_uniform_grid_points: " << 100.0*double(nb_point_with_noise)/cloud_point.size() << " % of points have been modified." << endl;

    cloud_point.swap(cloud_point_inter);
}

static void add_noise_to_points(const std::vector<Point3d>& candidates, std::vector<Point3d>& New_candidates)
{
    double x, y, z1, z2;
    unsigned int nbel = candidates.size();
    New_candidates.clear();
    //New_candidates.resize(nbel);
    srand((unsigned)time(NULL));
    for(unsigned int i=0; i<nbel; ++i)
    {
        gaussianNoise(x, y, false);
        gaussianNoise(z1, z2, false);
        Vector Noise(x,y,z1+z2);
        New_candidates.push_back(Point3d(candidates[i]+0.01*Noise));
        //assert(New_candidates[i]==Point3d(candidates[i]+0.01*Noise)); // ok
        //if(New_candidates[i]!=Point3d(candidates[i]+0.01*Noise))
        //{
        //    cout << "New_candidates[i].x() = " << New_candidates[i].x() << " New_candidates[i].y() = " << New_candidates[i].y() << " New_candidates[i].z() = " << New_candidates[i].z() << endl;
        //    cout << "Point3d(candidates[i]+0.01*Noise).x() = " << Point3d(candidates[i]+0.01*Noise).x() << " Point3d(candidates[i]+0.01*Noise).y() = " << Point3d(candidates[i]+0.01*Noise).y() << " Point3d(candidates[i]+0.01*Noise).z() = " << Point3d(candidates[i]+0.01*Noise).z() << endl;

        //    exit(-1);
        //}
    }
}

/*
    POINT SET READING (FROM FILE) AND SAVING (TO FILE)
*/
static void read_2Dpoints_from_file(const char* point_file,
                                    std::vector< Point3d >& cloud_point,
                                    bool add_noise=false)
{
    ifstream file(point_file, ios::in);  // on ouvre le fichier en lecture
    if(!file)
    {
        cerr << "Cannot open file: " << point_file << "!" << endl;
        //exit(-1);
        return;
    }

    double xmin=DBL_MAX, xmax=DBL_MIN, ymin=DBL_MAX, ymax=DBL_MIN;
    string line;
    while( getline(file, line) )
    { // for each point
        istringstream iss(line);
        double d1, d2;

        iss >> d1; iss >> d2;

        if(d1<xmin) xmin=d1;
        if(d1>xmax) xmax=d1;

        if(d2<ymin) ymin=d2;
        if(d2>ymax) ymax=d2;

        Point3d P(d1,d2,0);
        cloud_point.push_back(P);
    }
    // close the flux:
    file.close();

    if(add_noise)
        add_noise_to_points(cloud_point, xmin, xmax, ymin, ymax);
}

static void write_2Dpoints_to_file(const std::vector< Point3d >& cloud_point, const char* out_point_file)
{
    ofstream file(out_point_file, ios::out | ios::trunc);  // on ouvre le fichier en lecture
    if(!file)
    {
        cerr << "Cannot open file: " << out_point_file << "!" << endl;
        exit(-1);
    }

    std::vector< Point3d >::const_iterator its(cloud_point.begin()), ite(cloud_point.end());
    for(;its!=ite;++its)
    {
        file << its->x() << " " << its->y() << endl;
    }

    // close the flux:
    file.close();
}

#endif // Point_reading_generating_saving_H
