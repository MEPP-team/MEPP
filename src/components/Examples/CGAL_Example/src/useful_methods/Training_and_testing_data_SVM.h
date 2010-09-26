///////////////////////////////////////////////////////////////////////////
// Author: Vincent Vidal
// Year: 2010
// INSA-Lyon, LIRIS UMR 5205, M2DISCO.
///////////////////////////////////////////////////////////////////////////
#ifndef Training_and_testing_data_SVM_H
#define Training_and_testing_data_SVM_H

#include <string>

#include "../../../../../mepp/Polyhedron/polyhedron.h"

#include "../triangle_manager/CGAL_includes_and_types.h"

#include "../polyhedron_tools/Polyhedron_local_computation.h"

#include "../gnuplot/gnuplot_interface.h"

#if defined(linux) || defined(__linux) || defined(__linux__) || defined(__CYGWIN__)
	#ifndef UNIX_ARCHITECTURE
	#define UNIX_ARCHITECTURE
	#endif
#endif

static inline bool isFeatureEdgeBasedVertexColor(const Halfedge_handle pEdge,
                                                 char type=0)
{
    float   v1_r = pEdge->vertex()->color(0),
            v1_g = pEdge->vertex()->color(1),
            v1_b = pEdge->vertex()->color(2),
            v2_r = pEdge->prev()->vertex()->color(0),
            v2_g = pEdge->prev()->vertex()->color(1),
            v2_b = pEdge->prev()->vertex()->color(2);

    //cout << " v1_r = " << v1_r << " v1_g = " << v1_g << " v1_b = " << v1_b;
    //cout << " v2_r = " << v2_r << " v2_g = " << v2_g << " v2_b = " << v2_b << endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // we know that this type of edge cannot be feature edges!!
    if( cosDihedralAngle(pEdge)>.94 )
        return false;
    //////////////////////////////////////////////////////////////////////////////////////////////////
    if(!type)
        return (v1_r + v1_g + v1_b < 0.1) && (v2_r + v2_g + v2_b < 0.1);
    else
        return (v1_r>.4 || v1_g>.58) && (v2_r>.4 || v2_g>.58);
}
//////////////////////////////////////////////////////////////////////////////////////////
// recall (true detected true among all true)
static inline double recall(unsigned int nb_T_detected_T,
                            unsigned int nb_T_detected_F,
                            unsigned int nb_F_detected_T,
                            unsigned int nb_F_detected_F)
{
    return double(nb_T_detected_T)/std::max<double>(double(nb_T_detected_T+nb_T_detected_F),1.0);
}

// precision (true detected true among all detected true)
// it is also called the confidence [of the rule]
static inline double precision( unsigned int nb_T_detected_T,
                                unsigned int nb_T_detected_F,
                                unsigned int nb_F_detected_T,
                                unsigned int nb_F_detected_F)
{
    return double(nb_T_detected_T)/std::max<double>(double(nb_T_detected_T+nb_F_detected_T),1.0);
}

// Laplace correction smoothes probability estimates when the nb of instances
// covered by a rule is small
static inline double Laplace_corrected_confidence(  unsigned int nb_T_detected_T,
                                                    unsigned int nb_T_detected_F,
                                                    unsigned int nb_F_detected_T,
                                                    unsigned int nb_F_detected_F)
{
    return double(nb_T_detected_T+1.0)/double(nb_T_detected_T+nb_F_detected_T+2.0);
}

// Balanced Error Rate (BER)
static inline double BER(   unsigned int nb_T_detected_T,
                            unsigned int nb_T_detected_F,
                            unsigned int nb_F_detected_T,
                            unsigned int nb_F_detected_F)
{
    return 0.5*( double(nb_T_detected_F)/std::max<double>(double(nb_T_detected_F+nb_T_detected_T),1.0) + double(nb_F_detected_T)/std::max<double>(double(nb_F_detected_T+nb_F_detected_F),1.0) );
}

// for ROC curve
static inline double False_Positive_rate(   unsigned int nb_T_detected_T,
                                            unsigned int nb_T_detected_F,
                                            unsigned int nb_F_detected_T,
                                            unsigned int nb_F_detected_F)
{
    return double(nb_F_detected_T)/std::max<double>(double(nb_F_detected_T+nb_F_detected_F),1.0);
}

// it is the same as the recall
static inline double True_Positive_rate(   unsigned int nb_T_detected_T,
                                           unsigned int nb_T_detected_F,
                                           unsigned int nb_F_detected_T,
                                           unsigned int nb_F_detected_F)
{
    return double(nb_T_detected_T)/std::max<double>(double(nb_T_detected_T+nb_T_detected_F),1.0);
}

static inline double Error_rate(   unsigned int nb_T_detected_T,
                                   unsigned int nb_T_detected_F,
                                   unsigned int nb_F_detected_T,
                                   unsigned int nb_F_detected_F)
{
    return double(nb_F_detected_T+nb_T_detected_F)/std::max<double>(double(nb_T_detected_T+nb_T_detected_F+nb_F_detected_T+nb_F_detected_F),1.0);
}

static inline double Accuracy( unsigned int nb_T_detected_T,
                               unsigned int nb_T_detected_F,
                               unsigned int nb_F_detected_T,
                               unsigned int nb_F_detected_F
                             )
{
    return 1.0-Error_rate(nb_T_detected_T, nb_T_detected_F, nb_F_detected_T, nb_F_detected_F);
}
//////////////////////////////////////////////////////////////////////////////////////////
static inline double curvature( const Vertex_handle pVertex,
                                unsigned int which=2,
                                bool laplacian_averaging=false, // to use in case of noise
                                double radius = 0.07
                              )
{
    if(laplacian_averaging)
    {
        Halfedge_around_vertex_circulator h = pVertex->vertex_begin();
        Halfedge_around_vertex_circulator he = h;
        double sum=0.0;
        unsigned int cpt=0;
        CGAL_For_all(h,he)
        {
            sum+=curvature(h->prev()->vertex(),which);
            cpt++;
        }

        if(cpt) sum = sum/cpt;

        return sum;
    }
    else
    {
        switch(which)
        {
            case 0: // not discriminative
                return pVertex->closestObs()->Kmin; // detect corners
            case 1:
                return abs(pVertex->closestObs()->Kmin); // detect corners
            case 2: // can be discriminative
                return pVertex->closestObs()->Kmax; // detect lines of maximal curvatures
            case 3:
                return abs(pVertex->closestObs()->Kmax); // detect lines of maximal curvatures
            case 4: // can be discriminative
                return pVertex->closestObs()->Kmin * pVertex->closestObs()->Kmax; // well-known def of Gaussian curvature;detect corners
            case 5:
                return abs(pVertex->closestObs()->Kmin * pVertex->closestObs()->Kmax); // detect corners
            case 6: // can be discriminative
                return 0.5*(pVertex->closestObs()->Kmin + pVertex->closestObs()->Kmax); // well-known def of mean curvature; detect lines of maximal curvatures
            case 7:
                return abs(0.5*(pVertex->closestObs()->Kmin + pVertex->closestObs()->Kmax)); // detect lines of maximal curvatures
            case 8: // not very discriminative
                return pVertex->closestObs()->Kmax - pVertex->closestObs()->Kmin; // detect lines of maximal curvatures without corners
            case 9:
                return meanDihedralAngle(pVertex);
            case 10:
                return maxDihedralAngle(pVertex);
            case 11:
                return meanDirOneSidedCosTurningDirKmin(pVertex);
            case 12:
                return meanCurvatureAtVertex(pVertex);
            case 13:
                return G(pVertex, 1, radius);
            case 14:
                return multipleScaleSimilarity(pVertex, radius);
            default:
                return normal_variation(pVertex); // used in "Geometric snakes for triangular meshes"
        }
    }
}

static inline void curvature_measures(const Vertex_handle pVertex,
                                      vector< double >& row,
                                      bool laplacian_averaging=false)
{
    for(unsigned int i=0; i<11; i+=2)
        row.push_back(curvature(pVertex,i,laplacian_averaging));
}

static inline void curvature(   const Halfedge_handle pEdge,
                                vector< double >& row,
                                unsigned int which=2,
                                bool laplacian_averaging=false,
                                double radius = 0.07
                            )
{
    double tmpMean = 0.5*(curvature(pEdge->vertex(),which,laplacian_averaging)+curvature(pEdge->prev()->vertex(),which,laplacian_averaging));
    //double tmpO = abs(curvature(pEdge->vertex(),which,laplacian_averaging)-tmpMean);
    row.push_back(tmpMean);
    //row.push_back(tmpO);
}

static inline void curvature_measures(const Halfedge_handle pEdge,
                                        vector< double >& row,
                                        bool laplacian_averaging=false)
{
    for(unsigned int i=0; i<11; i+=2) // 6 times in loop
        curvature(pEdge, row, i, laplacian_averaging);
}

// Absolute curvature difference measure for 2 neighboring edges
// it is usefull to measure the local variation
static inline double curvature_diff_measures(   const Halfedge_handle pEdge,
                                                const Halfedge_handle pEdge2,
                                                vector< double >& row,
                                                bool laplacian_averaging=false)
{
    assert(pEdge!=Halfedge_handle());
    assert(pEdge2!=Halfedge_handle());
    for(unsigned int i=0; i<9; i+=2)
    {
        double tmpMean = 0.5*(curvature(pEdge->vertex(),i,laplacian_averaging)+curvature(pEdge->prev()->vertex(),i,laplacian_averaging));
        double tmpMean2 = 0.5*(curvature(pEdge2->vertex(),i,laplacian_averaging)+curvature(pEdge2->prev()->vertex(),i,laplacian_averaging));
        double tmpO = abs(curvature(pEdge->vertex(),i,laplacian_averaging)-tmpMean);
        double tmpO2 = abs(curvature(pEdge2->vertex(),i,laplacian_averaging)-tmpMean2);
        row.push_back(abs(tmpMean-tmpMean2)); row.push_back(abs(tmpO-tmpO2));
    }
}


// Training and testing data file generation
// The format of training and testing data file is:
// <label> <index1>:<value1> <index2>:<value2> ...
// <index> is an integer starting from 1 and <value> is a real number.
// Indices must be in ASCENDING order.
// Sample:
// +1 1:0.708333 2:1 3:1 4:-0.320755 5:-0.105023 6:-1 7:1 8:-0.419847 9:-1 10:-0.225806 12:1 13:-1
// -1 1:0.583333 2:-1 3:0.333333 4:-0.603774 5:1 6:-1 7:1 8:0.358779 9:-1 10:-0.483871 12:-1 13:1
// +1 1:0.166667 2:1 3:-0.333333 4:-0.433962 5:-0.383562 6:-1 7:-1 8:0.0687023 9:-1 10:-0.903226 11:-1 12:-1 13:1
// Each line contains an instance and is ended by a '\n' character.
// Labels in the testing file are only used to calculate accuracy or errors.
// If they are unknown, just fill the first column with any numbers.


static void fill_file_SVM(  ofstream& flux,
                            const std::vector< std::vector< double > >& data,
                            unsigned int max_row = 50000)
{
    unsigned int cpt=0;
    std::vector< vector< double > >::const_iterator it(data.begin()), ite(data.end());
    for(;(it!=ite)&&(cpt<max_row);++it)
    { // for each feature line:
        std::vector< double >::const_iterator it2(it->begin()), it2e(it->end());
        unsigned int index=1;
        if(*it2>=1) flux << "+";
        flux << int(*it2++); // the first row data represents the label!
        for(;it2!=it2e;++it2)
        {
            flux << " " << index << ":" << *it2;
            index++;
        }
        flux << endl;

        cpt++;
    }

}

static void generate_training_and_testing_files_SVM(const vector< vector< double > >& training_data,
                                                    const vector< vector< double > >& testing_data ,
                                                    bool concat_data )
{
    ofstream flux("./svm_training_data", ios::out | ((concat_data)?ios::app:ios::trunc));
    // open the file:
    if(!flux)
    {
        cout << "Cannot open the file: ./svm_training_data!!!" << endl;
        exit(-1);
    }

    fill_file_SVM(flux, training_data, 1000000);

    flux.close();
    //////////////////////////////////////////////////////////////////////
    ofstream flux2("./svm_testing_data", ios::out | ((concat_data)?ios::app:ios::trunc));
    // open the file:
    if(!flux)
    {
        cout << "Cannot open the file: ./svm_testing_data!!!" << endl;
        exit(-1);
    }

    fill_file_SVM(flux2, testing_data);

    flux2.close();

    cout << "generate_training_and_testing_files_SVM done!" << endl;
}

static void setFeatureEdgesFromVertexColor(PolyhedronPtr pMesh, char type=0)
{
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
        {

            if( isFeatureEdgeBasedVertexColor(pEdge,type) )
            {
                pEdge->SetFeature(true);
                pEdge->opposite()->SetFeature(true);
            }
            else
            {
                pEdge->SetFeature(false);
                pEdge->opposite()->SetFeature(false);
            }
        }
}

static void setFeatureEdgesFromVertexDirKmin(PolyhedronPtr pMesh)
{
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
        {

            if( 0.5*(oneCosTurningDirKmin(pEdge)+oneCosTurningDirKmin(pEdge->opposite())) > 0.991 )
            {
                pEdge->SetFeature(true);
                pEdge->opposite()->SetFeature(true);
            }
            else
            {
                pEdge->SetFeature(false);
                pEdge->opposite()->SetFeature(false);
            }
        }
}


static void setFeatureEdgesFromPrediction(PolyhedronPtr pMesh, const char* file_name="./svm_testing_data.predict") // 2 labels prediction only!!
{
    ifstream flux(file_name, ios::in);
    // open the file:
    if(!flux)
    {
        cout << "Cannot open the file: " << file_name << "!!!" << endl;
        exit(-1);
    }
    unsigned int cpt=0;
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
        {
            string line;

            if(!getline(flux, line)) break; // break is possible when there is not enough predicted data row

            istringstream iss(line);
            double d1;
            iss >> d1;

            if(d1>=1)
            {
                pEdge->SetFeature(true);
                pEdge->opposite()->SetFeature(true);
            }
            else
            {
                pEdge->SetFeature(false);
                pEdge->opposite()->SetFeature(false);
            }

            cpt++;
        }

    flux.close();

    if(!cpt) cout << "WARNING: setFeatureEdgesFromPrediction: svm_testing_data.predict file is empty!" << endl;
}

static inline void changeRhoModelFile(  double rho_init,
                                        double rho,
                                        const char* file_name="./svm_training_data.model")
{
    std::string Buffer(""); //Variable contenant le texte à réécrire dans le fichier
    std::vector <string> Vec_string(2500);
    std::ifstream ReadFile(file_name);
    if (ReadFile) //Si le fichier est trouvé
    {
        std::string line;
        int Line = 0;
        while ( std::getline( ReadFile, line ) ) //on parcours le fichier et on initialise line à la ligne actuelle
        {
            Line++;
            if(Line != 6)
            { //Si la ligne atteinte est différente de la ligne à supprimer...
                Buffer = line + "\n"; //On ajoute le contenu de la ligne dans le contenu à réécrire
            }
            else
            {
                std::ostringstream out;
                out << rho_init+rho;
                Buffer="rho ";
                Buffer+=out.str();
                Buffer+="\n";
            }

            Vec_string.push_back(Buffer);
            //cout << Vec_string.capacity() << endl;
        }
    }
    ReadFile.close(); //On ferme le fichier en lecture

    std::ofstream WriteFile( file_name ); //On ouvre ce même fichier en écriture
    std::vector <string>::const_iterator it(Vec_string.begin()), ite(Vec_string.end());
    for(;it!=ite; ++it)
        WriteFile << *it; //On écris le texte dedans

    WriteFile.close(); //et on ferme le fichier
}

// the output results won't be the same as the direct prediction given by the model
// because the scaling of the features is not necessary [-1,1]

static inline void setFeatureEdgesFromModel(PolyhedronPtr pMesh, bool scale_data=true) // 2 labels prediction only!!
{
#ifdef UNIX_ARCHITECTURE
    // we first generate the svm_testing_data.predict file using the SVM model
    if(scale_data)
    {
        FILE * ss = popen("svm-scale -l -1 -u 1 -s svm_training_data.range svm_training_data > svm_training_data.scale","r");
        if(ss == NULL)
        {
            fprintf(stderr, "Oops, I can't find %s.", "svm-scale");
            return;//exit(EXIT_FAILURE);
        }

        // terminer l'envoi de commande et fermer gnuplot
        pclose(ss);

        FILE * ss2 = popen("svm-scale -r svm_training_data.range svm_testing_data > svm_testing_data.scale","r");
        if(ss2 == NULL)
        {
            fprintf(stderr, "Oops, I can't find %s.", "svm-scale");
            return;//exit(EXIT_FAILURE);
        }

        // terminer l'envoi de commande et fermer gnuplot
        pclose(ss2);
    }

    system("svm-predict svm_testing_data.scale svm_training_data.model svm_testing_data.predict");

    // then we use the prediction results
    setFeatureEdgesFromPrediction(pMesh);
#else
	std::cout << "setFeatureEdgesFromModel does not work under Windows architecture." << std::endl;
#endif
}

/**
    this function generate one testing file for each mesh model (SVM data) and
    the training file corresponding to the other mesh features merging (SVM data)

    WORKS ONLY UNDER LINUX
**/
/*
static void generate_model_training_and_testing_files()
{
    ofstream flux

    FILE * ss = popen("svm-scale -l -1 -u 1 -s svm_training_data.range svm_training_data > svm_training_data.scale","r");
    if(ss == NULL)
    {
        fprintf(stderr, "Oops, I can't find %s.", "svm-scale");
        return;//exit(EXIT_FAILURE);
    }

    // terminer l'envoi de commande et fermer gnuplot
    pclose(ss);
}
*/

/**
    compareFeatureEdgesOnMeshAndGroundTruth
    :this method compares the current set of feature edges of pMesh
    with the one specified in a ground truth file; the nb of edges of
    pMesh must be equal to the nb of lines of file_ground_truth file;
    the result is expressed in terms of:
    - Banlanced Error Rate
    - recall
    - and precision

    The use_vertex_color flag permits to set the mesh edges
    to feature or not depending on its two vertices color
**/
static void compareFeatureEdgesOnMeshAndGroundTruth(PolyhedronPtr pMesh,
                                                    unsigned int& nb_T_T,
                                                    unsigned int& nb_T_F,
                                                    unsigned int& nb_F_F,
                                                    unsigned int& nb_F_T,
                                                    const char* file_ground_truth="./svm_training_data",
                                                    bool use_vertex_color=false
                                                    )
{
    ifstream flux(file_ground_truth, ios::in);//flux("./svm_training_data", ios::in);
    // open the file:
    if(!flux)
    {
        cout << "Cannot open the file: " << file_ground_truth << endl;
        exit(-1);
    }

    nb_T_T=nb_T_F=nb_F_F=nb_F_T=0;

    unsigned int cpt=0, nb_edges = CGAL::iterator_distance(pMesh->edges_begin(),pMesh->edges_end());

    if(use_vertex_color) // to set feature edges based on the vertex color
        setFeatureEdgesFromVertexColor(pMesh);

    // Here we are sure that the mesh feature edges are set
    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
        {
            string line;

            if(!getline(flux, line)) break; // break is possible when there is not enough predicted data row

            istringstream iss(line);
            double d1;
            iss >> d1;

            if(d1>=1)
            { // ground truth is true
                if(pEdge->isFeature())
                    nb_T_T++; // T_T
                else
                    nb_T_F++; // T_F
            }
            else
            { // ground truth is false
                if(pEdge->isFeature())
                    nb_F_T++; // F_T
                else
                    nb_F_F++; // F_F
            }

            cpt++;
        }

        if(cpt!=nb_edges)
            cout << "WARNING: compareFeatureEdgesOnMeshAndGroundTruth: not the same nb of lines and edges!!" << endl;
        else
        {
            cout << "Current recall (true detected true among all true) is " << recall(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;
            cout << "Current precision (true detected true among all detected true) is " << precision(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;
            cout << "Laplace corrected confidence is " << Laplace_corrected_confidence(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;

            cout << "False Positive Rate is " << False_Positive_rate(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;
            cout << "True Positive Rate is " << True_Positive_rate(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;

            cout << "Error Rate (ER) is " << Error_rate(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;
            cout << "Balanced Error Rate (BER) is " << BER(nb_T_T, nb_T_F, nb_F_T, nb_F_F) << endl;
        }
}


static void generate_training_and_testing_files_SVM(PolyhedronPtr pMesh,
                                                    PolyhedronPtr pObsMesh,
                                                    double min_dim,
                                                    bool use_vertex_color=true,
                                                    char type=0,
                                                    bool laplacian_averaging=false)
{
    vector< vector< double > > training_data, testing_data;
    bool isFirstLoop = true;

    for (   Edge_iterator pEdge = pMesh->edges_begin();
            pEdge != pMesh->edges_end();
            pEdge++)
    {
        vector< double > row;

        if( use_vertex_color && isFeatureEdgeBasedVertexColor(pEdge,type) ) // the edge is a sharp edge!!
        {
            row.push_back(1);
        }
        else
        {
            if(use_vertex_color || !pEdge->isFeature())
                row.push_back(-1);
            else
                row.push_back(1);
        }
        ////////////////////////////////////////////////////////////////////
        // geometric features:
        ////////////////////////////////////////////////////////////////////
        // dihedral angles
        row.push_back(cosDihedralAngle(pEdge)); // the more discriminative feature
        row.push_back(extendedCosDihedralAngle(pEdge)); // sensitive to noise!! the second more discriminative feature
        //row.push_back(regionBasedCosDihedralAngle(pEdge, 0.07*min_dim)); // less sensitive to noise

        curvature(pEdge,row,15); // normal variation
        ////////////////////////////////////////////////////////////////////
        // turning angles => degrades corner sharp edge detection and favor false-true when three edges are aligned
        //Halfedge_handle pE_out; // closest to pi
        //row.push_back(cosTurningEdgeAngleToTheStartorEnd(pEdge,true, pE_out));
        //curvature_diff_measures(pEdge, pE_out, row);
        //row.push_back(cosTurningEdgeAngleToTheStartorEnd(pEdge,false, pE_out));
        //curvature_diff_measures(pEdge, pE_out, row);
        ////////////////////////////////////////////////////////////////////
        //meanDirOneSidedCosTurningDirKmin(pVertex)
        ////////////////////////////////////////////////////////////////////
        training_data.push_back(row);
        testing_data.push_back(row); // used for prediction!
    }
    ///////////////////////////////////////////////////////////////////////
    double MinNrmMinCurvature,MaxNrmMinCurvature,MinNrmMaxCurvature,MaxNrmMaxCurvature;
    //0.04; in [0.01;0.4]
    // I have seen that a r greater than 0.24 was unuseful for SVM and current model set
    for(double r=0.01; r<=0.24; r+=0.04) // 6 times in loop
    {
        cout << "r = " << r << endl;

        principal_curvature(pObsMesh, true, r*min_dim, MinNrmMinCurvature, MaxNrmMinCurvature, MinNrmMaxCurvature, MaxNrmMaxCurvature);

        unsigned int cpt=0;
        for (   Edge_iterator pEdge = pMesh->edges_begin();
                pEdge != pMesh->edges_end();
                pEdge++)
        {
            ////////////////////////////////////////////////////////////////////
            // curvature tensor eigen values
            curvature_measures(pEdge, training_data[cpt], laplacian_averaging);
            curvature_measures(pEdge, testing_data[cpt], laplacian_averaging);
            ////////////////////////////////////////////////////////////////////
            //if(cpt==0) cout << training_data[cpt].size() << endl;
            ////////////////////////////////////////////////////////////////////
            training_data[cpt].push_back(regionBasedCosDihedralAngle(pEdge, r*min_dim)); // region dihedral angle
            testing_data[cpt].push_back(regionBasedCosDihedralAngle(pEdge, r*min_dim)); // region dihedral angle
            ////////////////////////////////////////////////////////////////////
            //if(cpt==0) cout << training_data[cpt].size() << endl;
            ////////////////////////////////////////////////////////////////////
            curvature(pEdge, training_data[cpt], 14, false, r*min_dim); // similarity
            curvature(pEdge, testing_data[cpt], 14, false, r*min_dim); // similarity
            ////////////////////////////////////////////////////////////////////
            //if(cpt==0) cout << training_data[cpt].size() << endl;
            ////////////////////////////////////////////////////////////////////
            cpt++;
        }
    }

    generate_training_and_testing_files_SVM(training_data,
                                            testing_data,
                                            false);
    // reset initial values:
    principal_curvature(pObsMesh, true, 0.01*min_dim, MinNrmMinCurvature, MaxNrmMinCurvature, MinNrmMaxCurvature, MaxNrmMaxCurvature);// to reset the initial curvatures!

}


static inline void generateROCcurvePNG( const vector < double > & False_Positive_values,
                                        const vector < double > & True_Positive_values,
                                        const char* out_file_name="./res_gnuplot/images/resRoc.png",
                                        bool remove_inter_data_file=true,
                                        const char* temp_file_folder="./res_gnuplot",
                                        const char* temp_file_code="temp")
{
    unsigned int nb_el = std::min<unsigned int>(False_Positive_values.size(),True_Positive_values.size());

    char current_file_name [256];

    sprintf (current_file_name, "%s/simulation_roc_curves_%s.txt",temp_file_folder,temp_file_code);

    ofstream flux(current_file_name, ios::out | ios::trunc);
    // open the file:
    if(!flux)
    {
        cout << "Cannot open the file: " << current_file_name << endl;
        exit(-1);
    }

    for(unsigned int i=0; i<nb_el; ++i)
    {
        if(False_Positive_values[i]>1.0)
            cout << "generateROCcurvePNG: WARNING False_Positive_values[i]=" << False_Positive_values[i] << endl;
        if(True_Positive_values[i]>1.0)
            cout << "generateROCcurvePNG: WARNING True_Positive_values[i]=" << True_Positive_values[i] << endl;

        if(False_Positive_values[i]<0.0)
            cout << "generateROCcurvePNG: WARNING negative False_Positive_values[i]=" << False_Positive_values[i] << endl;
        if(True_Positive_values[i]<0.0)
            cout << "generateROCcurvePNG: WARNING negative True_Positive_values[i]=" << True_Positive_values[i] << endl;

        flux << False_Positive_values[i] << " " << True_Positive_values[i];
        if(i<nb_el-1) flux << endl;
    }

    flux.close();

    gnuplot_interface g(current_file_name, "ROC", "False Positive Rate", "True Positive Rate","lines");
    g.gnuplot_png(out_file_name,-0.1,1.1, 0.0, 1.2);

    if(remove_inter_data_file)
    {
        char command [256];
        sprintf (command, "rm %s", current_file_name);
        system(command);
    }
}

static inline void readROCcurveFile(const char* file_name,
                                    vector < double > & False_Positive_values,
                                    vector < double > & True_Positive_values,
                                    bool extend_it = false // useful for models presenting at least one element of each class
                                    )
{
    False_Positive_values.clear(); False_Positive_values.reserve(32);
    True_Positive_values.clear(); True_Positive_values.reserve(32);

    ifstream flux(file_name, ios::in);
    // open the file:
    if(!flux)
    {
        cout << "readROCcurveFile: cannot open the file: " << file_name << "!!!" << endl;
        exit(-1);
    }

    bool first_line = true;
    string line;
    double d1, d2;
    while (getline(flux, line))
    {
        istringstream iss(line);
        iss >> d1; iss >> d2;

        assert(d1<=1.0);
        assert(d2<=1.0);

        assert(d1>=0.0);
        assert(d2>=0.0);

        if(first_line && extend_it && (d1<1.0 || d2<1.0))
        {
            if(d2!=0.0)
            {
                False_Positive_values.push_back(1.0);
                True_Positive_values.push_back(1.0);
            }
            else if(d1<1.0) // when there is not positive case at all (e.g. Sphere and torus)
            {
                False_Positive_values.push_back(1.0);
                True_Positive_values.push_back(0.0);
            }
        }

        False_Positive_values.push_back(d1);
        True_Positive_values.push_back(d2);

        first_line = false;
    }

    if(extend_it && (d1>0.0 || d2>0.0))
    {
        False_Positive_values.push_back(0.0);
        True_Positive_values.push_back(0.0);
    }
}

// the OX and OY data must be sorted (in increasing order or decreasing order, but consistently)
static inline double Area_Under_Curve(  const vector< double >& abs_X,
                                        const vector< double >& abs_Y)
{
    unsigned int nbe=abs_X.size();

    assert(nbe==abs_Y.size());

    if(nbe<=1) return 0.0;

    //cout << " nbe = " << nbe << endl;

    bool isIncresingOrder=true; // if false decreasing order

    unsigned k=0;
    while ( (k<nbe-1) && (abs_X[k]==abs_X[k+1]) && (abs_Y[k]==abs_Y[k+1]) ) k++;

    if(k<nbe-1) isIncresingOrder = !((abs_X[k]>abs_X[k+1]) || (abs_Y[k]>abs_Y[k+1]));

    double sum=0.0;
    if(isIncresingOrder)
    {
        //cout << " Increasing order..." << endl;
        for(unsigned int i=0; i<nbe-1; ++i)
        {
            //cout << "i = " << i << endl;
            //cout << "abs_X[i] = " << abs_X[i] << " abs_X[i+1] = " << abs_X[i+1] << endl;
            assert(abs_X[i+1]>=abs_X[i]);
            //cout << "abs_Y[i] = " << abs_Y[i] << " abs_Y[i+1] = " << abs_Y[i+1] << endl;
            assert(abs_Y[i+1]>=abs_Y[i]);
            sum += (abs_X[i+1]-abs_X[i]) * ( abs_Y[i] + 0.5*(abs_Y[i+1]-abs_Y[i]) );
        }
    }
    else
    {
        //cout << " Decreasing order..." << endl;
        for(unsigned int i=nbe-1; i>0; --i)
        {
            //cout << "i = " << i << endl;
            //cout << "abs_X[i] = " << abs_X[i] << " abs_X[i-1] = " << abs_X[i-1] << endl;
            assert(abs_X[i-1]>=abs_X[i]);
            //cout << "abs_Y[i] = " << abs_Y[i] << " abs_Y[i-1] = " << abs_Y[i-1] << endl;
            assert(abs_Y[i-1]>=abs_Y[i]);

            sum += (abs_X[i-1]-abs_X[i]) * ( abs_Y[i] + 0.5*(abs_Y[i-1]-abs_Y[i]) );
        }
    }

    return sum;
}

static inline double Area_Under_Curve(const char* file_name)
{
    vector < double > False_Positive_values, True_Positive_values;

    readROCcurveFile(file_name, False_Positive_values, True_Positive_values, true);

    return Area_Under_Curve(False_Positive_values, True_Positive_values);
}

/**
    generateMeanROCcurveFromIndividualROCcurves

 This method will read each individual ROC curve, will
 then extend them to start at (0,0) and to end at (1,1)

 Afterwards it will generate a png file corresponding to
 these extended ROC curves (and will also extend it).
 Curve files will be saved in the right directory (automatically).

 Vertical averaging:
 Then it will sample each extended ROC curve at regular interval
 and will generated a mean ROC curve taking the mean of the sampled values

     WORKS ONLY UNDER LINUX
**/
static inline void generateMeanROCcurveFromIndividualROCcurves( bool noisy=false,
                                                                bool minROCcurveInstead=false,
                                                                const char* common_name = "./res_gnuplot/data_simu_SVM_ROC_loc/simulation_roc_curves_whole_SVM_30_"
                                                                )
{
    unsigned int nb_samples = 400; // it seems to be sufficient
    double False_Positive_lb = 0.0, False_Positive_ub = 1.0;
    double False_Positive_step = (False_Positive_ub-False_Positive_lb)/double(nb_samples);

    const char* file_names[] = {  "1232_joint"
		                  ,"cone"
		                  ,"cup"
		                  ,"cut_cone"
		                  ,"cylinder"
		                  ,"fandisk"
		                  ,"renaud2"
		                  ,"PerfectSphere"
		                  ,"torus"
		                };

    unsigned int nb_files=7;
    ////////////////////////////////////////////////////////////
    string folder_name(common_name), category_name(common_name);
    int where_last_slash=1, i=folder_name.find('/');
    while(i!=-1)
    {
        where_last_slash=i;
        i=folder_name.find('/',i+1);
    }

    folder_name=folder_name.substr(0,where_last_slash);

    category_name=category_name.substr(category_name.find("simulation_roc_curves_whole_")+28);

    cout << "category = " << category_name.c_str() << endl;
    ////////////////////////////////////////////////////////////
    // a vector of local curves
    vector < vector < double > > False_Positive_values(nb_files);
    vector < vector < double > > True_Positive_values(nb_files);

    // for the mean ROC curve
    vector < double > False_Positive_values_G(nb_samples, 0.0);
    vector < double > True_Positive_values_G(nb_samples, 0.0);

    char current_file_name [256];
    cout << "generateMeanROCcurveFromIndividualROCcurves: nb_files = " << nb_files << endl;
    for(unsigned int i=0; i<nb_files; ++i)
    {
        if(noisy)
            sprintf (current_file_name, "%snoisy_%s.txt",common_name, file_names[i]);
        else
            sprintf (current_file_name, "%s%s.txt",common_name, file_names[i]);

        cout << "current_file_name = " << current_file_name;
        readROCcurveFile(   current_file_name,
                            False_Positive_values[i],
                            True_Positive_values[i],
                            true);

        cout << " AUC = " << Area_Under_Curve(False_Positive_values[i], True_Positive_values[i]) << endl;

        if(noisy)
            sprintf (current_file_name, "%s/resRoc_noisy_%s.png",folder_name.c_str(),file_names[i]);
        else
            sprintf (current_file_name, "%s/resRoc_%s.png",folder_name.c_str(),file_names[i]);

        generateROCcurvePNG(False_Positive_values[i],
                            True_Positive_values[i],
                            current_file_name);
    }
    cout << "generateMeanROCcurveFromIndividualROCcurves: individual curves have been generated." << endl;
    // for each false positive value, we will compute the mean of all ROC curves
    for(double f_Pos=1.0; f_Pos>=0.0; f_Pos-=False_Positive_step)
    {
        unsigned int index_l, index_r; // closest index values of false positive

        unsigned int index_g = int(f_Pos/False_Positive_step);
        assert(index_g<nb_samples);
        //cout << "global current index = " << index_g << endl;
        False_Positive_values_G[index_g]=f_Pos;

        //cout << "False_Positive_values_G[index_g] = " << False_Positive_values_G[index_g] << endl;

        if(minROCcurveInstead)
            True_Positive_values_G[index_g]=1.0;

        for(unsigned int i=0; i<nb_files; ++i)
        {
            unsigned int nb_max = False_Positive_values[i].size();
            index_l = index_r = 0; // remind that values are stored in reverse order
            //cout << "init index_l = " << index_l << endl;

            // we are looking for the fisrt index such as
            // False_Positive_values[i][index_l]<=f_Pos
            while(  (index_l<nb_max) &&
                    (False_Positive_values[i][index_l]>f_Pos) // for a fixed fp rate, we need
                                                              // to take the maximum tp rate
                                                              // with that fp rate
                 )
            {
                //cout << "False_Positive_values_G[index_g] = " << False_Positive_values_G[index_g] << " False_Positive_values[i][index_l] = " << False_Positive_values[i][index_l] << endl;
                index_l++;
            }

            //cout << "index_l = " << index_l << endl;
            double current_true_pos_val;
            if( (index_l<nb_max) &&
                (False_Positive_values[i][index_l]==f_Pos)
              )
            {
                current_true_pos_val = False_Positive_values[i][index_l];
            }
            else if( (index_r < index_l) && (index_l<nb_max))
            {
                index_r = index_l;
                index_l--;
                double diff = abs(False_Positive_values[i][index_r]-False_Positive_values[i][index_l]);

                current_true_pos_val =  abs(False_Positive_values[i][index_r]-f_Pos)/diff*True_Positive_values[i][index_l] // when f_Pos==False_Positive_values[i][index_l] we must have False_Positive_values[i][index_l] at the end!
                                        +
                                        abs(False_Positive_values[i][index_l]-f_Pos)/diff*True_Positive_values[i][index_r];

            }
            else
            {
                //cout << "boundary case";
                if(index_l>=nb_max)
                { // case where the local ROC curve is strictly greater than f_Pos at the left
                  // => interpolate between (0,0) and the left-most point of local ROC curve
                    //cout << " at the left";
                    index_r = index_l = nb_max-1;
                    double diff = False_Positive_values[i][index_l];

                    current_true_pos_val = f_Pos/diff * True_Positive_values[i][index_l];

                    //cout << "diff = " << diff << " current_true_pos_val = " << current_true_pos_val << endl;
                }
                else
                { // case where the local ROC curve is less than f_Pos at the right
                  // => interpolate between (1,1) and the right-most point of local ROC curve
                    //cout << " at the right";
                    index_r = index_l;
                    double diff = 1.0-False_Positive_values[i][index_l];
                    current_true_pos_val =  abs(1.0-f_Pos)/diff*True_Positive_values[i][index_l]
                                            +
                                            abs(False_Positive_values[i][index_l]-f_Pos)/diff;
                }

                //cout << endl;
            }

            if(current_true_pos_val>1.0)
                cout << "current_true_pos_val = " << current_true_pos_val << endl;

            assert(current_true_pos_val<=1.0);

            if(minROCcurveInstead)
                True_Positive_values_G[index_g] = std::min<double>(current_true_pos_val, True_Positive_values_G[index_g]);
            else
                True_Positive_values_G[index_g] = (i*True_Positive_values_G[index_g]+current_true_pos_val)/double(i+1.0);
            //cout << "True_Positive_values_G[index_g] = " << True_Positive_values_G[index_g] << endl;
        }
    }

    //
    if((False_Positive_values_G[0]>0.0) || (True_Positive_values_G[0]>0.0))
    {
        False_Positive_values_G.insert(False_Positive_values_G.begin(),0.0);
        True_Positive_values_G.insert(True_Positive_values_G.begin(),0.0);

        nb_samples++;
    }

    if((False_Positive_values_G[nb_samples-1]<1.0) || (True_Positive_values_G[nb_samples-1]<1.0))
    {
        False_Positive_values_G.push_back(1.0);
        True_Positive_values_G.push_back(1.0);

        nb_samples++;
    }

    cout << "MeanRoc AUC = " << Area_Under_Curve(False_Positive_values_G, True_Positive_values_G) << endl;

    char extra_suffix [256];
    if(noisy)
    {
        sprintf (current_file_name, "%s/resNoisyMeanRoc.png",folder_name.c_str());

        sprintf (extra_suffix, "%snoisy_glob",category_name.c_str());

        generateROCcurvePNG(False_Positive_values_G,
                            True_Positive_values_G,
                            current_file_name,
                            false,
                            "./res_gnuplot/Mean_ROC_curve_glob",
                            extra_suffix);
    }
    else
    {
        sprintf (current_file_name, "%s/resMeanRoc.png",folder_name.c_str());

        sprintf (extra_suffix, "%sinit_glob",category_name.c_str());

        generateROCcurvePNG(False_Positive_values_G,
                            True_Positive_values_G,
                            current_file_name,
                            false,
                            "./res_gnuplot/Mean_ROC_curve_glob",
                            extra_suffix);
    }

    cout << "generateMeanROCcurveFromIndividualROCcurves: global curve has been generated." << endl;
}
#endif // Training_and_testing_data_SVM_H
