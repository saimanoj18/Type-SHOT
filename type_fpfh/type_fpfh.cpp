#include<type_rops/typerops_bits.h>

/* callback for qsort(): */
typedef struct {float delta; int i;} D;

// comparing function
int D_cmp (const void *a_, const void *b_)
{
    D* a = (D *)a_;
    D* b = (D *)b_;
    return (a->delta > (b->delta)? 1: ((a->delta < b->delta)? -1: 0));
}

/* limits: */
#define M_MAX   100    /* max number of dimensions */
#define N_MAX   100   /* max value of common denominator */

/* Pascal's triangle: */
static int *binomial[N_MAX+1], b_data[(N_MAX+1) * (N_MAX+2) / 2];

/*
 * Initialize quantization module.
 */
void t_init ()
{
    int n, k, *b = b_data;
    for (n=0; n<=N_MAX; n++) {
        binomial[n] = b; b += n + 1;            /* row pointer */
        binomial[n][0] = binomial[n][n] = 1;    /* 1st & last coeffs */
        for (k=1; k<n; k++)                     /* compute coeffs in the middle */
            binomial[n][k] = binomial[n-1][k-1] + binomial[n-1][k];
    }
}


/*
 * Compute number of reconstruction points.
 * Input:
 *   m - number of dimensions
 *   n - precision parameter
 * Returns:
 *   number of types with m dimensions and denominator n
 */
int t_points (int m, int n)
{
    assert(m < M_MAX);
    assert(n < N_MAX);
    assert(n+m-1 < N_MAX);
    assert(binomial[0][0] == 1);
    /* return number of types: */
    return binomial[n+m-1][m-1];
}


int t_quant (int m, int n, std::vector<float>& p)
{
    int i, j, idx, n1, Delta, s;
    int k[M_MAX]; D d[M_MAX];
    //cout << "p :" << p[0] << endl;
    /* quantize distribution p[] to type k[]: */
    for (n1=0,i=0; i<m; i++) n1 += (k[i] = (int) floor(p[i] * n + 0.5));     // unconstrained quantization
    Delta = n1 - n;
    if (Delta != 0) {
        for (i=0; i<m; i++) {d[i].delta = (float)k[i] - p[i] * n; d[i].i = i;} // compute & sort errors
        std::qsort(d, m, sizeof(D), D_cmp);
        if (Delta > 0) {for (j=m-Delta; j<m; j++) k[d[j].i] --;}               // ensure that sum_i k[i] = n
        else           {for (j=0; j<abs(Delta); j++) k[d[j].i] ++;}
    }

    /* compute type index: */
    for (i=0, idx=0; i<m-2; i++) {
        s = 0;
        for (j=0; j<k[i]; j++)
            s += binomial[n-j+m-i-2][m-i-2];
        idx += s;
        n -= k[i];
    }
    idx += k[m-2];

    return idx;
}

/*
 * Reconstruct quantized probability distribution.
 * Input:
 *   m - number of dimensions, n - precision parameter
 *   idx - index of a reconstruction point
 * Output:
 *   p[0..m-1] - reconstructed probability distribution
 * Returns:
 *   1 - success; 0 - invalid index
 */
int t_reconst (int m, int n, int idx, std::vector<float>& p)
{
    p.resize(m);
    int i, j, k[M_MAX], s, x; float n_inv = 1.f/(float)n;

    /* check index: */
    if (idx < 0 || idx >= t_points(m,n))
        return 0;

    /* decode type: */
    for (i=0; i<m-2; i++)
    {
        s = 0;
        for (j=0; j<n; j++)
        {
            x = binomial[n-j+m-i-2][m-i-2];
            if (idx - s < x) break;
            s += x;
        }
        k[i] = j;
        idx -= s;
        n -= j;
    }
    k[m-2] = idx;
    k[m-1] = n - idx;

    //cout<< "idx : " << idx << endl;
    /* convert it to distribution: */
    for (j=0; j<m; j++)
    {
        p[j] = (float)k[j] * n_inv;
        //cout << p[j] << " ";
    }

    return 1;
}

void create_LUT (int m, int n, Eigen::MatrixXf &LUT)
{
    int range = t_points(m,n);
    cout << "\nNo. of Lattices : " <<range << endl;
    LUT.resize(range, range);
    std::vector<float> p;
    std::vector<float> q;



    for (int idx_lut_i = 0; idx_lut_i < range; idx_lut_i++)
    {
        t_reconst(m, n, idx_lut_i, p);

        for (int idx_lut_j = 0; idx_lut_j < range; idx_lut_j++)
        {
            t_reconst(m, n, idx_lut_j, q);

            float dist = 0;
            for (int dist_idx = 0; dist_idx < m; dist_idx++)
            {
                dist = dist + sqrt((p[dist_idx] - q[dist_idx]) * (p[dist_idx] - q[dist_idx]));
                //dist = dist + fabs((p[dist_idx] - q[dist_idx]));

            }
            LUT(idx_lut_i, idx_lut_j) = dist;
            //cout << "i : " << idx_lut_i << " j : " << idx_lut_j << "dist : " << dist << endl;
        }

    }
}

void normalize(std::vector<float>& V, std::vector<float>& out)
{
    int range = V.size();
    out.resize(range);
    float sum = 0;
    float neg = 0;
    float check = 0;
    for (int i = 0; i < range; i++)
        if (V[i] <= neg)
        {
            neg = V[i];
            check = 1;
        }
    if (check == 1)
        for (int i = 0 ; i < range; i++)
            V[i] = V[i] +abs(neg);


    for (int i = 0; i < range; i++)
        sum = sum + V[i];
    for (int i = 0; i < out.size(); i++)
        if(sum > 0)
            out[i] = (float)V[i]/(float)sum;

}


/*************************************************
 *
 *
 Implementation of Lattice Quantizer based on the steps mentioned in PhD
 thesis of Vijay Chandrashekar titled,
 "Low-bitrate Image Retrieval With Compressed Histogram of Gradients Descriptors"

 m, n and P are the input parameters

 m represents the number of bins in the input histogram
 n is the quatization parameter
 P is the normalized histogram
 *
 *
 * ************************************************/




void compute_compressed_fpfh_from_FPFH(pcl::PointCloud<pcl::FPFHSignature33 >& shot_descriptors_here, std::vector<desc>& CBSHOT_descriptors, int m, int n, int grids, int bins)
{
    for (int i = 0; i < (int)shot_descriptors_here.size(); i++)
    {
        desc descriptor; descriptor.values.resize(grids);

        for (int j = 0 ; j < grids ; j++)
        {
            std::vector<float> temp, out;

            for (int k = 0; k < bins; k++)
            {
                float x;
                x = shot_descriptors_here[i].histogram[ j*bins + k ];

                if(pcl_isnan(x))
                {
                    temp.push_back(0);
                    cout << x << endl;
                }
                else if (pcl_isnan(-x))
                    temp.push_back(0);
                else temp.push_back(x);
            }

            normalize(temp, out);

            float output = 0; Eigen::VectorXf k;
            output = t_quant(m, n, out);
            descriptor.values[j] = output;
        }
        CBSHOT_descriptors.push_back(descriptor);
    }

}


int main()
{

    int m, n, grids, bins;
    m = 5;
    n = 8;
    grids = 7;
    bins = 5;

    cbshot cb;
    // Read a PCD file from disk.

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/scene030_0.pcd", cb.cloud2);
    //pcl::io::loadPolygonFilePLY("../sample_files/scene030_0.ply", cb.mesh2);

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/Squirrel003_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/Squirrel003_0.ply", cb.mesh1);

    pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/scene005_0.pcd", cb.cloud2);
    //pcl::io::loadPolygonFilePLY("../sample_files/scene005_0.ply", cb.mesh2);

    pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/Doll018_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/Doll018_0.ply", cb.mesh1);

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/mario000_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/mario000_0.ply", cb.mesh1);

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/PeterRabbit001_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/PeterRabbit001_0.ply", cb.mesh1);

    t_init();
    Eigen::MatrixXf LUT;
    cout << "Started to construct Look Up Table" << endl;
    create_LUT(m, n, LUT);
    cout << "LUT constructed! " << endl;
/*
    //cout << LUT << endl;
    ofstream myfile;
      myfile.open ("example.txt");
      myfile << LUT;
      myfile.close();
*/
    cb.calculate_normals (0.02);

    //    Eigen::Matrix4f temp;
    //    temp << 0.995665, 	-0.00287336,	-0.0929723,	-0.0966733,
    //            0.00810565,	0.998401,	0.0559452,	-0.0522716,
    //            0.0926635,	-0.0564565,	0.994095,	-0.0271486,
    //            0,	0	,0,	1;

    //    cb.calculate_voxel_grid_keypoints_for_evaluation(0.01, temp);

    cb.calculate_voxel_grid_keypoints (0.01);

    cb.get_keypoint_indices();

    //calculate_my_don2_keypoints (cloud1, cloud1_normals, cloud2, cloud2_normals, cloud1_keypoints, cloud2_keypoints);

    //calculate_my_super_duper_keypoints (cloud1, cloud1_normals, cloud2, cloud2_normals, cloud1_keypoints, cloud2_keypoints);

    //cb.calculate_SHOT(0.10);
    //cb.calculate_rops(0.04);
    cb.calculate_FPFH(0.05);


    compute_compressed_fpfh_from_FPFH( cb.cloud1_fpfh, cb.compressed_shot1, m, n, grids, bins);
    compute_compressed_fpfh_from_FPFH( cb.cloud2_fpfh, cb.compressed_shot2, m, n, grids, bins);


    /**************************************************/


    clock_t start2, end2;
    double cpu_time_used2;
    start2 = clock();


    pcl::Correspondences corrs;

    Eigen::MatrixXf dist;
    dist.resize(cb.compressed_shot1.size(), cb.compressed_shot2.size());
    dist.setZero();

    for (int i = 0; i < cb.compressed_shot1.size(); i++)
    {
        for (int j = 0; j < cb.compressed_shot2.size(); j++)
        {
            for (int l = 0; l < grids; l++)
                dist(i,j) = dist(i,j) + LUT(cb.compressed_shot1[i].values[l], cb.compressed_shot2[j].values[l]);
            //dist(i,j) = dist(i,j) + (((cb.compressed_shot1[i].values[l]- cb.compressed_shot2[j].values[l]))*((cb.compressed_shot1[i].values[l]- cb.compressed_shot2[j].values[l])));
        }
    }




    for (int i = 0; i < dist.rows(); i++)
    {
        float small = dist(i,0);
        int index_i = i;
        int index_j = 0;
        for(int j = 0; j < dist.cols(); j++)
        {
            if (dist(i,j) < small)
            {
                small = dist(i,j);
                index_i = i;
                index_j = j;
            }
        }
        int temp = 0;
        for (int j = 0; j < dist.rows(); j++)
        {
            if (dist(j,index_j) < small)
                temp = 1;
        }
        if (temp == 0)
        {
            pcl::Correspondence corr;
            //cout << "here " << endl;
            //cout << cb.cloud1_keypoints_indices[index_i] << endl;
            corr.index_query = cb.cloud1_keypoints_indices[index_i];// vulnerable
            corr.index_match = cb.cloud2_keypoints_indices[index_j];// vulnerable
            corr.distance = dist(index_i,index_j);

            corrs.push_back(corr);
        }

    }



    end2 = clock();
    cpu_time_used2 = ((double) (end2 - start2)) / CLOCKS_PER_SEC;
    std::cout << "Time taken for Type-SHOT NN Matching : " << (double)cpu_time_used2 << std::endl;

    /************************************************************/



    cout << "No. of Reciprocal Correspondences : " << corrs.size() << endl;



    pcl::CorrespondencesConstPtr corrs_const_ptr = boost::make_shared< pcl::Correspondences >(corrs);

    pcl::Correspondences corr_shot;
    pcl::registration::CorrespondenceRejectorSampleConsensus< pcl::PointXYZ > Ransac_based_Rejection_shot;
    Ransac_based_Rejection_shot.setInputSource(cb.cloud1.makeShared());
    Ransac_based_Rejection_shot.setInputTarget(cb.cloud2.makeShared());
    Ransac_based_Rejection_shot.setInlierThreshold(0.02);
    Ransac_based_Rejection_shot.setInputCorrespondences(corrs_const_ptr);
    Ransac_based_Rejection_shot.getCorrespondences(corr_shot);

    cout << "Mat : \n" << Ransac_based_Rejection_shot.getBestTransformation()<< endl;

    cout << "True correspondences after RANSAC : " << corr_shot.size() << endl;


    boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer (new pcl::visualization::PCLVisualizer ("3D Viewer"));
    viewer->setBackgroundColor (255, 255, 255);



    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color1(cb.cloud1.makeShared(), 255, 0, 0);
    viewer->addPointCloud<pcl::PointXYZ> (cb.cloud1.makeShared(), single_color1, "sample cloud1");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "sample cloud1");
    //viewer->addCoordinateSystem (1.0);
    viewer->initCameraParameters ();

    Eigen::Matrix4f t;
    t<<1,0,0,0.6,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1;

    //cloudNext is my target cloud
    pcl::transformPointCloud(cb.cloud2,cb.cloud2,t);

    //int v2(1);
    //viewer->createViewPort (0.5,0.0,0.1,1.0,1);
    pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZ> single_color2(cb.cloud2.makeShared(), 0, 0, 255);
    viewer->addPointCloud<pcl::PointXYZ> (cb.cloud2.makeShared(), single_color2, "sample cloud2");
    viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 4, "sample cloud2");



    viewer->addCorrespondences<pcl::PointXYZ>(cb.cloud1.makeShared(), cb.cloud2.makeShared(), /*corrs*/ corr_shot, "correspondences"/*,v1*/);



    while (!viewer->wasStopped ())
    {
        viewer->spinOnce (100);
        boost::this_thread::sleep (boost::posix_time::microseconds (100000));
    }




    return 0;
}

