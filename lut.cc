#include "pegasus/pegasus_runtime.h"
#include "pegasus/timer.h"

using namespace std;
using namespace seal;
using namespace gemini;
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>

void load_matrix(const std::string& filename, std::vector<std::vector<double>>& matrix) {
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        int numRows = 0;
        int numCols = 0;
        
        // 统计行数和列数
        while (std::getline(file, line)) {
            numRows++;
	    numCols = 0;
            std::istringstream iss(line);
            double value;
            while (iss >> value) {
                numCols++;
            }
        }
        printf("numRows:%d, numCols:%d\n",numRows,numCols);
        // 重新定位文件指针到文件开头
        file.clear();
        file.seekg(0, std::ios::beg);
        
        // 创建矩阵
        matrix.resize(numRows, std::vector<double>(numCols, 0));
        
        // 读取数据到矩阵
        int row = 0, col = 0;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            double value;
            while (iss >> value) {
                matrix[row][col] = value;
                col++;
            }
            row++;
            col = 0;
        }
        
        file.close();
    }
}


void matrixMultiply(std::vector<std::vector<double>>& result, const std::vector<std::vector<double>>& mat1, const std::vector<std::vector<double>>& mat2) {
    // 获取矩阵的维度
    int n = mat1.size();
    int m = mat1[0].size();
    int p = mat2[0].size();

    // 初始化结果矩阵为0


    // 进行矩阵乘法
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            for (int k = 0; k < m; ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

}

void encode_and_encrypte_expand(PegasusRunTime& pg_rt, vector<vector<Ciphertext>> & encrypted_expand, const std::vector<std::vector<double>> plain_expand, size_t num_vertex, size_t num_packed)
{
    auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < num_vertex; i++){
	vector<Plaintext> plain_matrix1(num_vertex/num_packed);
	vector<Ciphertext> encrypted_matrix1(num_vertex);
	for (size_t j = 0; j < num_vertex/num_packed; j++)
	{
		auto start = std::chrono::high_resolution_clock::now();
		vector<double> pack_v;//(poly_modulus_degree/2);
		for (size_t k = 0; k < num_packed; k++)
		{
			size_t index_now = i*num_vertex + j*num_packed +k;
			pack_v.insert(pack_v.end(), plain_expand[index_now].begin(), plain_expand[index_now].end());
		}
   /*
		encoder.encode(pack_v, scale, plain_matrix1[j]);
		encryptor.encrypt(plain_matrix1[j], encrypted_matrix1[j]);*/
   		pg_rt.EncodeThenEncrypt(pack_v, encrypted_matrix1[j]);
   

	}
	encrypted_expand [i] = encrypted_matrix1;
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "total encode and encrypte time: " << elapsed.count() << " s\n";

}

void Aggregation(PegasusRunTime& pg_rt, vector<Ciphertext>& encrypted_product, const vector<vector<Ciphertext>>& encrypted_adj_expand, const vector<vector<Ciphertext>>& encrypted_node_expand, size_t num_vertex,size_t num_packed)
{
    auto start = std::chrono::high_resolution_clock::now();	
    for (size_t i = 0; i < num_vertex; i++)
    {
        vector<Ciphertext> row_products(num_vertex/num_packed);
        for (size_t j = 0; j < num_vertex/num_packed; j++)
        {
	    auto start = std::chrono::high_resolution_clock::now();
            /*
            evaluator.multiply(encrypted_adj_expand[i][j], encrypted_node_expand[i][j], row_products[j]);
            evaluator.relinearize_inplace(row_products[j], relin_keys);*/
	   Ciphertext* add_row_products = & row_products[j];
            pg_rt.runtime_->MulRelin(add_row_products,encrypted_adj_expand[i][j],encrypted_node_expand[i][j]);
	    pg_rt.runtime_->RescaleNext(add_row_products);
            //evaluator.rescale_to_next_inplace(row_products[j]);
	        // 记录结束时间点
    	    auto finish = std::chrono::high_resolution_clock::now();
    // 计算代码执行的时间
            std::chrono::duration<double> elapsed = finish - start;
    	    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
	    std::cout << "total time: " << elapsed.count()*num_vertex*num_vertex/num_packed << " s\n";
       

        }
	if(i == 0)
	  for (size_t j = 0; j < num_vertex/num_packed; j++)
          {
            encrypted_product[j] = row_products[j];
          }
        else
          for (size_t j = 0; j < num_vertex/num_packed; j++)
          {
            pg_rt.Add(encrypted_product[j], row_products[j]);
          }
    }
    auto finish = std::chrono::high_resolution_clock::now();
    // 计算代码执行的时间
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "total Aggregation time: " << elapsed.count() << " s\n";

}

void ReLU(PegasusRunTime& pg_rt, PegasusRunTime::Parms pp, Ciphertext ckks_ct, std::vector<double> &slots)
{
    double s2c_time{0.};
    AutoTimer timer1(&s2c_time);
    CHECK_AND_ABORT(pg_rt.SlotsToCoeffs(ckks_ct));
    timer1.stop();
    printf("S2c time without extract %f sec\n",s2c_time/1000);
    std::vector<lwe::Ctx_st> lwe_ct;
    double s2c_ks_time{0.};
    {
      AutoTimer timer(&s2c_ks_time);
      CHECK_AND_ABORT(pg_rt.ExtraAllCoefficients(ckks_ct, lwe_ct));
    }

    double accum_s2c_err{0.};
    for (int i = 0; i < pp.nslots; ++i) {
      accum_s2c_err += std::abs(slots[i] - pg_rt.DecryptLWE(lwe_ct[i]));
    }
    printf("\nS2C.error 2^%f\t%f sec\n", std::log2(accum_s2c_err / pp.nslots), s2c_ks_time / 1000.);

    std::string tag;
    std::function<double(double)> target_func;

    double lut_time{0.};

    AutoTimer timer(&lut_time);

    pg_rt.ReLU(lwe_ct.data(), lwe_ct.size());
    tag = "ReLU";
    target_func = [](double e) -> double { return std::max(0., e); };
    timer.stop();

    double accum_lut_err{0.};
    for (int i = 0; i < pp.nslots; ++i) {
      double gnd = target_func(slots[i]);
      double cmp = pg_rt.DecryptLWE(lwe_ct[i]);
      accum_lut_err += std::abs(gnd - cmp);
    }

     printf("LUT (%s).error 2^%f\t%f sec\n", tag.c_str(), std::log2(accum_lut_err / pp.nslots), lut_time / 1000.);
}

int main()
{

    
    
    using namespace gemini;
    PegasusRunTime::Parms pp;

    pp.lvl0_lattice_dim = lwe::params::n();
    pp.lvl1_lattice_dim = 1 << 12;
    pp.lvl2_lattice_dim = 1 << 16;
    pp.nlevels = 4; // CKKS levels
    pp.scale = std::pow(2., 40);
    pp.nslots = 1 << 10;
    pp.s2c_multiplier = 1.;

    PegasusRunTime pg_rt(pp, /*num_threads*/ 1);
    
    // 从txt文件读取edge数据
    std::string filename_edge_expand = "/home/wyt/projects/plain_gcn/Cora_edge_expand_256.txt";
    std::vector<std::vector<double>> adj_matrix;
    load_matrix(filename_edge_expand, adj_matrix);

    
    // 从txt文件读取node数据
    std::string filename_node = "/home/wyt/projects/plain_gcn/Cora_node_f16_ex_256_expand.txt";
    std::vector<std::vector<double>> node_matrix;
    load_matrix(filename_node, node_matrix);    





    int num_vertex_sq = adj_matrix.size();
    int m1 = adj_matrix[0].size();
    printf("num_vertex_sq x m1  = %d x %d\n",num_vertex_sq,m1);
    int n2 = node_matrix.size();
    int f_len = node_matrix[0].size();
    printf("n2 x f_len  = %d x %d\n",n2,f_len);


   // compute in plaintext
    std::string filename_node_o = "/home/wyt/projects/plain_gcn/Cora_node_f16_ex_256.txt";
    std::string filename_edge = "/home/wyt/projects/plain_gcn/Cora_edge_256.txt";
    std::vector<std::vector<double>> node_matrix_o;
    load_matrix(filename_node_o, node_matrix_o);  

    size_t num_vertex = node_matrix_o.size();
    printf("num_vertex = %d \n",num_vertex);

    std::vector<std::vector<double>> adj_matrix_o;
    load_matrix(filename_edge, adj_matrix_o);
    std::vector<std::vector<double>> matrix_result_plaintext(num_vertex, std::vector<double>(f_len, 0));
    matrixMultiply(matrix_result_plaintext, adj_matrix_o, node_matrix_o);


    
 //encode and encrypte expand adj 
    size_t num_packed =  pp.nslots/2/f_len;  
    vector<vector<Ciphertext>> encrypted_adj_expand(num_vertex, vector<Ciphertext>(num_vertex/num_packed));
    encode_and_encrypte_expand(pg_rt, encrypted_adj_expand, adj_matrix,num_vertex, num_packed);
     printf("finish encode and encrypte adj\n");

 //encode and encrypte expand nodes 
    vector<vector<Ciphertext>> encrypted_node_expand(num_vertex, vector<Ciphertext>(num_vertex/num_packed));
    encode_and_encrypte_expand(pg_rt, encrypted_node_expand, node_matrix, num_vertex, num_packed);
    printf("finish encode and encrypte vertices\n");
 // Compute Combination
   //
 // Compute Aggregation
    vector<Ciphertext> encrypted_product(num_vertex/num_packed); //
    Aggregation(pg_rt, encrypted_product, encrypted_adj_expand, encrypted_node_expand, num_vertex, num_packed);
    printf("finish Aggregation\n");
    // encrypted_product is not in expand form
 // Compute ReLU



 // 计算并验证relu结果
 
    vector<Plaintext> plain_product(num_vertex/num_packed);
    vector<vector<double>> matrix_result_packed(num_vertex/num_packed);
    for (size_t i = 0; i < num_vertex/num_packed; i++)
    {
	pg_rt.DecryptThenDecode(encrypted_product[i],matrix_result_packed[i]);
	ReLU(pg_rt,pp,encrypted_product[i],matrix_result_packed[i]);
        //decryptor.decrypt(encrypted_product[i], plain_product[i]);
	//encoder.decode(plain_product[i], matrix_result_packed[i]);
    }
   printf("start check\n");
   printf("matrix_result_packed size : %d x %d\n",matrix_result_packed.size(),matrix_result_packed[0].size());

// 验证Aggregation结果
    double err = 0;
    for (size_t i = 0; i < num_vertex/num_packed; i++)
    {
	for (size_t j = 0; j < num_packed*f_len; j++){
	size_t squence = i*num_packed*f_len+ j;
        size_t row = squence/f_len;
        size_t column = squence%f_len;
        //printf("row x column: %d %d\n",row,column);
	 err += std::abs(matrix_result_plaintext[row][column] - matrix_result_packed[i][j]);
       //printf("%f : %f\n",matrix_result_plaintext[row][column],matrix_result_packed[i][j]);
	}
       
    }

   printf("err: %f\n",err);






    return 0;
}






