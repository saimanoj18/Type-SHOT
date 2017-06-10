#include<type_shot/typeshot_bits.h>


/* callback for qsort(): */
typedef struct {float delta; int i;} D;

//int D_cmp (const void *D l, const void *D r)
//{
//	return (l->delta > r->delta) ? 1 : ((l->delta < r->delta)? -1 : 0);
//}

int D_cmp (const void *a_, const void *b_)
{
D* a = (D *)a_;
D* b = (D *)b_;
  return (a->delta > (b->delta)? 1: ((a->delta < b->delta)? -1: 0));
}

/* limits: */
#define M_MAX   50    /* max number of dimensions */
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
    for (i=0; i<m-2; i++) {
        s = 0;
        for (j=0; j<n; j++) {
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

    /* convert it to distribution: */
    for (j=0; j<m; j++)
        p[j] = ((float)k[j]) * n_inv;

    return 1;
}



void create_LUT (int m, int n_, Eigen::MatrixXf &LUT)
{
    cout << t_points(m, n_) << endl;
    LUT.resize(t_points(m,n_), t_points(m,n_));
    std::vector<float> p;
    std::vector<float> q;


    int range = t_points(m,n_);
    for (int idx_lut_i = 0; idx_lut_i < range; idx_lut_i++)
    {
        t_reconst(m, n_, idx_lut_i, p);

        for (int idx_lut_j = 0; idx_lut_j < range; idx_lut_j++)
        {
            t_reconst(m, n_, idx_lut_j, q);

            float dist = 0;
            for (int dist_idx = 0; dist_idx < m; dist_idx++)
            {
                dist = dist + sqrt((p[dist_idx] - q[dist_idx]) * (p[dist_idx] - q[dist_idx]));
            }
            LUT(idx_lut_i, idx_lut_j) = dist;
            //cout << "i : " << idx_lut_i << " j : " << idx_lut_j << "dist : " << dist << endl;
        }

    }
}


void normalize(std::vector<float>& V, std::vector<float>& out)
{
    out.resize(V.size());
    float sum = 0;
    float neg = 0;
    float check = 0;
    for (int i = 0; i < V.size(); i++)
        if (V[i] <= neg)
        {
            neg = V[i];
            check = 1;
        }
    if (check == 1)
        for (int i = 0 ; i < V.size(); i++)
            V[i] = V[i] +abs(neg);


    for (int i = 0; i < V.size(); i++)
        sum = sum + V[i];
    for (int i = 0; i < out.size(); i++)
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
void compressedSHOT (int m, int n, std::vector<float>& moments, Eigen::VectorXf& k, float& output)
{
    //normalize the moment vector
    k.resize(moments.size());

    if(m < n)
    {
        //std::cout << "Error ... m is less than n" <<std::endl;
    }
    Eigen::VectorXf moment_vector;
    moment_vector.resize(moments.size());

    for (int i = 0; i < moments.size(); ++i)
    {
        moment_vector(i) = moments[i];
    }

    /*
    %%   step 1 from the thesis
    for i = 1 : m
        k_dash(i) = floor (( n*P(i) ) + 0.5);
    end

    n_dash = sum(k_dash);
    done = 0;
    */

    Eigen::VectorXf k_dash,delta; k_dash.resize(m);k_dash.resize(m);delta.resize(m);
    for (int i = 0; i < m; i++)
    {
        k_dash[i] = floor(( n * moment_vector(i) ) + 0.5);
    }

    float n_dash = k_dash.sum();
    float done = 0;

    /*
    %% Step 2 from the thesis
    if (n_dash == n)
        k = k_dash;
        done = 1;
    end

    if (done == 0)
        if (n_dash ~= n)
            for i = 1 : m
                delta(i) = k_dash(i) - ( n*P(i) );
            end
        end


        j_i(1:m,1:2) = 1;%initializing
        temp = 1;% temp
        for i = 1 : m
            j_i(i,temp) = [delta(i)];
            j_i(i,temp+1) = i;
        end

        j_i;

        sorted_j_i = sortrows(j_i, 1);

        */

    if (n_dash == n)
    {
        k = k_dash;
        done = 1;
    }
    if (done == 0)
    {
        if (n_dash != n)
        {
            for (int i = 0; i < m; i++)
            {
                delta(i) = k_dash(i) - ( n*moment_vector(i) );
            }
        }

        std::vector<mypair> vector_of_pairs;
        vector_of_pairs.resize(m);

        for (int k = 0 ; k < m ; k++)
        {
            vector_of_pairs[k].first = delta(k);
            vector_of_pairs[k].second = k;
        }

        std::sort (vector_of_pairs.begin(), vector_of_pairs.end(), comparator);

        /*

        %% Step 3----Mistake in d > 0 and d < 0 conditions...
        d = n_dash - n;

        if ( d > 0 )
            k = k_dash;
            for i = 1 : d
                k(sorted_j_i(m + 1 - i, 2)) = k(sorted_j_i(m + 1 - i, 2)) - 1;
            end
        end

        if ( d < 0 )
            k = k_dash;
            for i = 1 : abs(d)
                k(sorted_j_i( i, 2)) = k(sorted_j_i( i, 2)) + 1;
            end
        end
    end
        */

        float d = n_dash - n;

        if ( d > 0)
        {
            k = k_dash;
            for ( int i = 0; i < d; i++)
            {
                k(vector_of_pairs[m - 1 - i].second)--;
            }
        }

        if ( d < 0)
        {
            k = k_dash;
            for ( int i = 0; i < std::abs(d); i++)
            {
                k(vector_of_pairs[i].second) = k(vector_of_pairs[i].second) + 1;
            }
        }
    }

    /*
    %% Step 4 -- ToDO Enumeration of the Lattice points

    enumeration = 1;
    output = 0;
    for j = 1 : n-2
        for i = 0 : (k(j) - 1)
            k_l = 0;
            for l = 1 : (j - 1)
                k_l = k_l + k(l);
            end
            output = output + (factorial( m - j + (n - i - k_l) - 1)/( (factorial(n - i - k_l)) * (factorial(m-j-1)) ));
        end
    end
*/

    /* compute type index: */
    int s;
    for (int i=0, idx=0; i<m-2; i++)
    {
        s = 0;
        for (int j=0; j<k[i]; j++)
        {
            s += binomial[n-j+m-i-2][m-i-2];
        }
        idx += s;
        n -= k[i];
    }
    output += k[m-2];

    /*

    output = 0;
    for (int j = 1; j <= n-2; j++)
    {
        for ( int i = 0; i <= (k(j-1) -1); i++)
        {
            float k_l = 0;
            {
                for ( int l = 1; l <= (j - 1); l++)
                {
                    k_l = k_l + k(l-1);
                }
            }
            output = output + (factorial( m - j + (n - i - k_l) - 1)/( (factorial(n - i - k_l)) * (factorial(m-j-1)) ));
        }
    }
    */

}



void compute_compressed_shot_from_SHOT(pcl::PointCloud<pcl::SHOT352>& shot_descriptors_here, std::vector<desc>& CBSHOT_descriptors)
{
    //CBSHOT_descriptors.resize(shot_descriptors_here.size());
//0.12  0.1 0.05 0.01 0.09 0.32 0.11 0.2

    std::vector<float> temp, out;
    temp.resize(8);temp[0] = 12; temp[1] = 10; temp[2] = 5; temp[3] = 1; temp[4] =  9; temp[5]= 32; temp[6] = 11; temp[7] = 20;
    //temp.resize(5);temp[0] = 12; temp[1] = 28; temp[2] = 17; temp[3] = 27; temp[4] =  16;// temp[5]= 32; temp[6] = 11; temp[7] = 20;

//[12, 28, 17, 27, 16]
            normalize(temp, out);
            for (int idd = 0; idd < 8; idd++ )
                cout << out[idd] << endl;

            //for (int l = 0; l < 4; l++ )
            //cout << out[l] << endl;
            //cout << endl;
            float output = 0; Eigen::VectorXf k;
            output = t_quant(8, 8, out);

            cout << output << endl;
            //cout << endl;

            std::vector<float> q;
            t_reconst(8,8, output, q);

            for (int idd = 0; idd < 8; idd++ )
                cout << q[idd] << endl;

//[0.125, 0.125, 0.0625, 0.0625, 0.125, 0.1875, 0.125, 0.1875].


}



int main()
{
    cbshot cb;
    // Read a PCD file from disk.

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/scene030_0.pcd", cb.cloud2);
    //pcl::io::loadPolygonFilePLY("../sample_files/scene030_0.ply", cb.mesh2);

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/Squirrel003_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/Squirrel003_0.ply", cb.mesh1);

    pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/scene005_0.pcd", cb.cloud2);
    //pcl::io::loadPolygonFilePLY("../sample_files/scene005_0.ply", cb.mesh2);

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/Doll018_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/Doll018_0.ply", cb.mesh1);

    //pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/mario000_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/mario000_0.ply", cb.mesh1);

    pcl::io::loadPCDFile<pcl::PointXYZ>("../sample_files/PeterRabbit001_0.pcd", cb.cloud1);
    //pcl::io::loadPolygonFilePLY("../sample_files/PeterRabbit001_0.ply", cb.mesh1);

    t_init();
    Eigen::MatrixXf LUT;
    cout << "Started to construct Look Up Table" << endl;
    create_LUT(2,2, LUT);
    cout << "LUT constructed! " << endl;

    //cout << LUT << endl;

    cb.calculate_normals (0.02);//0.02(0.001UWA)

    cb.calculate_voxel_grid_keypoints (0.01);//0.01(0.005UWA)
//    Eigen::Matrix4f temp;
//    temp << 0.995665, 	-0.00287336,	-0.0929723,	-0.0966733,
//            0.00810565,	0.998401,	0.0559452,	-0.0522716,
//            0.0926635,	-0.0564565,	0.994095,	-0.0271486,
//            0,	0	,0,	1;

//    cb.calculate_voxel_grid_keypoints_for_evaluation(0.01, temp);

    cb.get_keypoint_indices();

    //calculate_my_don2_keypoints (cloud1, cloud1_normals, cloud2, cloud2_normals, cloud1_keypoints, cloud2_keypoints);

    //calculate_my_super_duper_keypoints (cloud1, cloud1_normals, cloud2, cloud2_normals, cloud1_keypoints, cloud2_keypoints);

    cb.calculate_SHOT(0.10);//0.10(0.03UWA)


//    for (int i = 200; i < 250; i++)
//    {
//        cout << endl;
//        for (int j = 0; j < 352; j++)
//        cout << cb.cloud1_shot[i].descriptor[j] << " ";
//    }

    compute_compressed_shot_from_SHOT( cb.cloud1_shot, cb.compressed_shot1);
    //compute_compressed_shot_from_SHOT( cb.cloud2_shot, cb.compressed_shot2);


    /**************************************************/

    pcl::Correspondences corrs;

    Eigen::MatrixXf dist;
    dist.resize(cb.cloud1_shot.size(), cb.cloud2_shot.size());
    dist.setZero();

    //    for (int i = 0; i < cb.histograms1.size(); i++)
    //    {
    //        for (int j = 0; j < cb.histograms2.size(); j++)
    //        {
    //            //dist(i,j) = (sh1.signature_descriptors[i] - sh2.signature_descriptors[j]).norm();

    //            Eigen::VectorXf one, two;
    //            one.resize(27);
    //            two.resize(27);

    //            for(int k = 0; k < 27; k++)
    //            {
    //                one[k] = cb.histograms1[i].histogram[k];
    //            }

    //            for(int k = 0; k < 27; k++)
    //            {
    //                two[k] = (cb.histograms2[j].histogram[k]);
    //            }

    //            dist(i,j) = (one-two).norm();


    //        }
    //    }


    for (int i = 0; i < cb.cloud1_shot.size(); i++)
    {
        for (int j = 0; j < cb.cloud2_shot.size(); j++)
        {
            //dist(i,j) = (sh1.signature_descriptors[i] - sh2.signature_descriptors[j]).norm();

            Eigen::VectorXf one, two;
            one.resize(352);
            two.resize(352);

            for(int k = 0; k < 352; k++)
            {
                one[k] = cb.cloud1_shot[i].descriptor[k];
            }

            for(int k = 0; k < 352; k++)
            {
                two[k] = cb.cloud2_shot[j].descriptor[k];
            }

            dist(i,j) = (one-two).norm();


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
        //int temp = 0;
        //for (int j = 0; j < dist.rows(); j++)
        //{
        //    if (dist(j,index_j) < small)
        //        temp = 1;
        //}
        //if (temp == 0)
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

