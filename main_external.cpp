#include "seal/seal.h"
using namespace std;
using namespace seal;
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

void encode_and_encrypte_expand(Encryptor& encryptor,CKKSEncoder& encoder, double scale, vector<vector<Ciphertext>> & encrypted_expand, const std::vector<std::vector<double>> plain_expand, size_t num_vertex, size_t num_packed)
{
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
		encoder.encode(pack_v, scale, plain_matrix1[j]);
		encryptor.encrypt(plain_matrix1[j], encrypted_matrix1[j]);
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "total time: " << elapsed.count()*num_vertex*num_vertex/num_packed << " s\n";
	}
	encrypted_expand [i] = encrypted_matrix1;
    }

}

void Aggregation(Evaluator& evaluator, RelinKeys& relin_keys, vector<Ciphertext>& encrypted_product, const vector<vector<Ciphertext>>& encrypted_adj_expand, const vector<vector<Ciphertext>>& encrypted_node_expand, size_t num_vertex,size_t num_packed)
{
    for (size_t i = 0; i < num_vertex; i++)
    {
        vector<Ciphertext> row_products(num_vertex/num_packed);
        for (size_t j = 0; j < num_vertex/num_packed; j++)
        {
	    auto start = std::chrono::high_resolution_clock::now();	
            evaluator.multiply(encrypted_adj_expand[i][j], encrypted_node_expand[i][j], row_products[j]);
            evaluator.relinearize_inplace(row_products[j], relin_keys);
            evaluator.rescale_to_next_inplace(row_products[j]);
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
            evaluator.add_inplace(encrypted_product[j], row_products[j]);
          }
    }

}
int main()
{
    // 从txt文件读取edge数据
    std::string filename_edge_expand = "Cora_edge_expand_1024.txt";
    std::vector<std::vector<double>> adj_matrix;
    load_matrix(filename_edge_expand, adj_matrix);

    
    // 从txt文件读取node数据
    std::string filename_node = "Cora_node_f16_ex_1024_expand.txt";
    std::vector<std::vector<double>> node_matrix;
    load_matrix(filename_node, node_matrix);    




    // 创建加密参数
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = 8192;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 40, 60 }));

    // 创建SEAL上下文
    auto context = SEALContext(parms);

    // 创建密钥
    KeyGenerator keygen(context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    SecretKey secret_key = keygen.secret_key();
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);

    // 创建加密器和解密器
    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);

    // 创建编码器和求值器
    CKKSEncoder encoder(context);
    Evaluator evaluator(context);
    double scale = pow(2.0, 40);

    
    

    int num_vertex_sq = adj_matrix.size();
    int m1 = adj_matrix[0].size();
    printf("num_vertex_sq x m1  = %d x %d\n",num_vertex_sq,m1);
    int n2 = node_matrix.size();
    int f_len = node_matrix[0].size();
    printf("n2 x f_len  = %d x %d\n",n2,f_len);


   // compute in plaintext for 
    std::string filename_node_o = "Cora_node_f16_ex_1024.txt";
    std::string filename_edge = "Cora_edge_1024.txt";
    std::vector<std::vector<double>> node_matrix_o;
    load_matrix(filename_node_o, node_matrix_o);  

    size_t num_vertex = node_matrix_o.size();
    printf("num_vertex = %d \n",num_vertex);

    std::vector<std::vector<double>> adj_matrix_o;
    load_matrix(filename_edge, adj_matrix_o);
    std::vector<std::vector<double>> matrix_result_plaintext(num_vertex, std::vector<double>(f_len, 0));
    matrixMultiply(matrix_result_plaintext, adj_matrix_o, node_matrix_o);


    
 //encode and encrypte expand adj 
    size_t num_packed =  poly_modulus_degree/2/f_len;  
    vector<vector<Ciphertext>> encrypted_adj_expand(num_vertex, vector<Ciphertext>(num_vertex/num_packed));
    encode_and_encrypte_expand(encryptor, encoder, scale, encrypted_adj_expand, adj_matrix,num_vertex, num_packed);
     printf("finish encode and encrypte adj\n");

 //encode and encrypte expand nodes 
    vector<vector<Ciphertext>> encrypted_node_expand(num_vertex, vector<Ciphertext>(num_vertex/num_packed));
    encode_and_encrypte_expand(encryptor, encoder, scale, encrypted_node_expand, node_matrix, num_vertex, num_packed);
    printf("finish encode and encrypte vertices\n");
 // Compute Combination
   /* not yet*/

 // Compute Aggregation
    vector<Ciphertext> encrypted_product(num_vertex/num_packed); //
    Aggregation(evaluator, relin_keys, encrypted_product, encrypted_adj_expand, encrypted_node_expand, num_vertex, num_packed);
    printf("finish Aggregation\n");
    // encrypted_product is not in expand form
 // Compute ReLU
 // 解密结果
    vector<Plaintext> plain_product(num_vertex/num_packed);
    vector<vector<double>> matrix_result_packed(num_vertex/num_packed);
    for (size_t i = 0; i < num_vertex/num_packed; i++)
    {
        decryptor.decrypt(encrypted_product[i], plain_product[i]);
	encoder.decode(plain_product[i], matrix_result_packed[i]);
    }
   printf("start check\n");
   printf("matrix_result_packed size : %d x %d\n",matrix_result_packed.size(),matrix_result_packed[0].size());

// 验证结果
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

