#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include<tuple>

float Determinant_3x3(const std::vector<std::vector<float>>& mat) {
    return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
        mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
        mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
}

// 4x4 matrisin determinantýný hesaplayan fonksiyon
float Determinant_4x4(const std::vector<std::vector<float>>& mat) {
    float det = 0.0f;

    // Ýlk satýr üzerinden kofaktör geniþlemesi
    for (int col = 0; col < 4; col++) {
        // 3x3 alt matrisi oluþtur
        std::vector<std::vector<float>> minor_matrix(3, std::vector<float>(3));
        for (int i = 1; i < 4; i++) {
            int minor_col = 0;
            for (int j = 0; j < 4; j++) {
                if (j == col) continue;
                minor_matrix[i - 1][minor_col] = mat[i][j];
                minor_col++;
            }
        }

        // Determinantý hesapla
        float minor_det = Determinant_3x3(minor_matrix);
        det += (col % 2 == 0 ? 1 : -1) * mat[0][col] * minor_det;
    }
    return det;
}

//Matrisleri toplamak için kullanýlýr
std::vector<std::vector<float>> Matrix_Sum(const std::vector<std::vector<float>>& Mat_1, const std::vector<std::vector<float>>& Mat_2) {
    // Boyut uyumsuzluðu kontrolü
    if (Mat_1.size() != Mat_2.size() || Mat_1[0].size() != Mat_2[0].size()) {
        throw std::invalid_argument("Matris boyutlarý uyumsuz.");
    }

    int row = Mat_1.size();
    int col = Mat_1[0].size();

    std::vector<std::vector<float>> result(row, std::vector<float>(col, 0));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            result[i][j] = Mat_1[i][j] + Mat_2[i][j];
        }
    }
    return result;
}

std::vector<std::vector<float>> Mat_4x4_Ters(const std::vector<std::vector<float>>& Mat) {
    int n = 4;
    if (Determinant_4x4(Mat) != 0) {
        //[A|I] matrisi oluþturulur.
        std::vector<std::vector<float>> Augmented_Matrix(n, std::vector<float>(2 * n, 0));

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Augmented_Matrix[i][j] = Mat[i][j];
            }
            Augmented_Matrix[i][n + i] = 1;
        }
        //Gauss-Jordan Eleme Yöntemi
        for (int i = 0; i < n; ++i) {
            float diagElement = Augmented_Matrix[i][i];
            if (diagElement == 0) {
                // Satýr deðiþimi uygulayýn
                for (int k = i + 1; k < n; ++k) {
                    if (Augmented_Matrix[k][i] != 0) {
                        std::swap(Augmented_Matrix[i], Augmented_Matrix[k]);
                        break;
                    }
                }
                diagElement = Augmented_Matrix[i][i];
            }
            for (int j = 0; j < 2 * n; ++j) {
                Augmented_Matrix[i][j] /= diagElement;
            }
            for (int k = 0; k < n; ++k) {
                if (k != i) {
                    float factor = Augmented_Matrix[k][i];
                    for (int j = 0; j < n * 2; ++j) {
                        Augmented_Matrix[k][j] -= factor * Augmented_Matrix[i][j];
                    }
                }
            }
        }
        //Ters matrisi oluþturulur.
        std::vector<std::vector<float>> Inverse_Matrix(n, std::vector<float>(n, 0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Inverse_Matrix[i][j] = Augmented_Matrix[i][n + j];
            }
        }

        return Inverse_Matrix;
    }
    else {
        std::cout << "Determinant sifir." << std::endl;
        return Mat;
    }
}

//Matrisin bir sayý ile çarpýlmasý için kullanýlýr.
std::vector<std::vector<float>> Multiply(const std::vector<std::vector<float>>& Mat, float a) {
    std::vector<std::vector<float>> result(Mat.size(), std::vector<float>(Mat[0].size(), 0));
    for (int i = 0; i < Mat.size(); i++) {
        for (int j = 0; j < Mat[0].size(); j++) {
            result[i][j] = Mat[i][j] * a;
        }
    }
    return result;
}

std::vector<std::vector<float>> CholeskyDecomposition(const std::vector<std::vector<float>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<float>> L(n, std::vector<float>(n, 0.0f));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            float sum = 0.0f;
            if (j == i) {
                for (int k = 0; k < j; k++)
                    sum += L[i][k] * L[i][k];
                if (matrix[i][i] - sum < 0) {
                    std::cout << "Error: negative value in sqrt: " << matrix[i][i] - sum << " i= " << i << std::endl;
                    return L;
                }
                L[i][i] = std::sqrtf(matrix[i][i] - sum);
            }
            else {
                for (int k = 0; k < j; k++)
                    sum += L[i][k] * L[j][k];
                if (L[j][j] == 0) {
                    std::cout << "Error: division by 0, L[" << j << "][" << j << "] = " << L[j][j] << " j= " << j << std::endl;
                    return L;
                }

                L[i][j] = (matrix[i][j] - sum) / L[j][j];
            }
        }
    }
    return L;
}

//Ýki matrisin çarpýmý için kullanýlýr.
std::vector<std::vector<float>> Matrix_Multiply(const std::vector<std::vector<float>>& Mat_1, const std::vector<std::vector<float>>& Mat_2) {
    int row = Mat_1.size();
    int col = Mat_2[0].size();
    int mid = Mat_2.size();

    std::vector<std::vector<float>> result(row, std::vector<float>(col, 0));
    if (Mat_1[0].size() != Mat_2.size()) {
        std::cout << "Matrislerin boyutlari uyumsuzdur." << std::endl;
        return result;
    }
    else {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                for (int k = 0; k < mid; k++) {
                    result[i][j] += Mat_1[i][k] * Mat_2[k][j];
                }
            }
        }
        return result;
    }
}

//Bir matrisin transpozunu almak için kullanýlýr.
std::vector<std::vector<float>> Transpose(const std::vector<std::vector<float>>& Mat) {
    std::vector<std::vector<float>> Mat_Tra(Mat[0].size(), std::vector<float>(Mat.size(), 0));
    for (int i = 0; i < Mat.size(); i++) {
        for (int j = 0; j < Mat[0].size(); j++) {
            Mat_Tra[j][i] = Mat[i][j];
        }
    }
    return Mat_Tra;
}

//4x1 -> 4x9 matris
std::vector<std::vector<float>> Genisletme(const std::vector<std::vector<float>>& Mat) {
    std::vector<std::vector<float>> result(4, std::vector<float>(9, 0));
    for (int i = 0; i < Mat.size(); i++) {
        for (int j = 0; j < 9; j++) {
            result[i][j] = Mat[i][0];
        }
    }
    return result;
}

//4x9 -> 4x1 matris
std::vector<std::vector<float>> Daraltma(const std::vector<std::vector<float>>& Mat) {
    std::vector<std::vector<float>> result(4, std::vector<float>(1, 0));
    for (int i = 0; i < Mat.size(); i++) {
        for (int j = 0; j < 9; j++) {
            result[i][0] = Mat[i][j];
        }
    }
    return result;
}

//Sigma noktalarýnýn hesaplanmasý(4 x 9 matris oluþturur)
std::vector<std::vector<float>> Sigma_Points(const std::vector<std::vector<float>>& Mat, const std::vector<std::vector<float>>& P) {
    int n = Mat.size();
    const float a = 1e-3, k = 0;
    float lambda = (a * a * (n + k)) - n;

    // Sigma noktalarý matrisini oluþtur (boyut: n x (2n+1))
    std::vector<std::vector<float>> sigma_points(n, std::vector<float>(2 * n + 1, 0));

    // Ýlk sigma noktasý = Ortalama durum vektörü
    for (int i = 0; i < n; i++) {
        sigma_points[i][0] = Mat[i][0];
    }

    // Kare kök matris hesaplama (Cholesky ayrýþýmý kullanýlmalý)
    std::vector<std::vector<float>> sqrt_matrix = CholeskyDecomposition(Multiply(P, (n + lambda))); // 4 x 4 matris

    // Sigma noktalarýný hesapla
    for (int i = 0; i < n; i++) { 
        for (int j = 0; j < n; j++) { 
            sigma_points[j][i + 1] = Mat[j][0] + sqrt_matrix[j][i] * sqrtf(n + lambda);  // Pozitif sigma noktasý
            sigma_points[j][i + 1 + n] = Mat[j][0] - sqrt_matrix[j][i] * sqrtf(n + lambda);  // Negatif sigma noktasý
           // std::cout << sqrt_matrix[i][j] << "     ";
        }
        //std::cout << std::endl;
    }
    return sigma_points;
}

//Sigma Noktalarýnýn Aðýrlýklarýnýn Hesaplanmasý
std::vector<std::vector<float>> Agýrlýk_Hesaplama(int L){
	const float a = 0.9, k = 0;
	float lambda = (a * a * (L + k)) - L;
	std::vector<std::vector<float>> W_c(2 * L + 1, std::vector<float>(1, 0));
	std::vector<std::vector<float>> W_m(2 * L + 1, std::vector<float>(1, 0));
	std::vector<std::vector<float>> W(2 * L + 1, std::vector<float>(2, 0));
    W_m[0][0] = lambda / (L + lambda);
    W_c[0][0] = lambda / (L + lambda) + (1 - (a * a) + k);
	for (int i = 1; i <= 2 * L; i++) {
		W_m[i][0] = 1 / (2 * (L + lambda));
        W_c[i][0] = W_m[i][0];
	}
    for (int i = 0; i <= 2 * L; i++) {
        W[i][0] = W_m[i][0];
        W[i][1] = W_c[i][0];
    }
	return W;
}

//Sigma Bazlý Kalman Filtresi Tahmin
std::tuple<std::vector<std::vector<float>>, std::vector<std::vector<float>>> Sigma_Kalman_Filter_Prediction(
    float delta_t, float x, float v_x, float y, float v_y, float a_x, float a_y) {
    ///Zaman adýmý
    float sigma_x2 = 10, sigma_xv2 = 10, cov_xv = 10;//rastgele belirlenmiþtir. 
    float sigma_y2 = 10, sigma_yv2 = 10, cov_yv = 10;

    float sigma_x2_k1{}, sigma_v2x_k1{}, cov_xv_k1{};//X 
    float sigma_y2_k1{}, sigma_v2y_k1{}, cov_yv_k1{};//Y 

	//Durum vektörünü ayarla
    std::vector<std::vector<float>> x_k = { {x},{v_x},{y},{v_y} };

    //Hata Kovaryans Matrisi(k zamanýndaki)
    std::vector<std::vector<float>> p_k = { {sigma_x2,cov_xv,0,0},{cov_xv,sigma_xv2,0,0},{0,0,sigma_y2,cov_yv},{0,0,cov_yv,sigma_yv2} };

    //Sigma noktalarýný ayarla
    std::vector<std::vector<float>> sigma_p = Sigma_Points(x_k, p_k); // 4 x 9 matris

    //Durum geçiþ fonksiyonunu ayarla
    std::vector<std::vector<float>> sigma_x_k1(sigma_p.size(), std::vector<float>(sigma_p[0].size(), 0)); // 4 x 9 matris

    for (int i = 0; i < sigma_p[0].size(); i++) {
        sigma_x_k1[0][i] = sigma_p[0][i] + sigma_p[1][i] * delta_t + 0.5f * a_x * delta_t * delta_t; //x1 = x0 + v_x0 * t + 0.5 * a_x * t^2
        sigma_x_k1[1][i] = sigma_p[1][i] + a_x * delta_t; //v_x1 = v_x0 + a_x * t

        sigma_x_k1[2][i] = sigma_p[2][i] + sigma_p[3][i] * delta_t + 0.5f * a_y * delta_t * delta_t; //y1 = y0 + v_y0 * t + 0.5 * a_y * t^2
        sigma_x_k1[3][i] = sigma_p[3][i] + a_y * delta_t; //v_y1 = v_y0 + a_y * t
    }

    //Durum tahmini (k + 1 zamanýndaki)
    std::vector<std::vector<float>> sigma_weights = Agýrlýk_Hesaplama(sigma_p.size());   // 9 x 2 matris
    std::vector<std::vector<float>> x_k1(x_k.size(), std::vector<float>(x_k[0].size(), 0)); // 4 x 1 matris
	std::vector<std::vector<float>> p_k1(p_k.size(), std::vector<float>(p_k[0].size(), 0)); // 4 x 4 matris

    for (int i = 0; i < sigma_weights.size(); i++) {
		x_k1[0][0] += sigma_weights[i][0] * sigma_x_k1[0][i];
		x_k1[1][0] += sigma_weights[i][0] * sigma_x_k1[1][i];
		x_k1[2][0] += sigma_weights[i][0] * sigma_x_k1[2][i];
		x_k1[3][0] += sigma_weights[i][0] * sigma_x_k1[3][i];
    }

    //Hata kovaryansý tahmini (k + 1 zamanýndaki)
    x_k1 = Genisletme(x_k1); //4x1 -> 4x9

    std::vector<std::vector<float>> diff = Matrix_Sum(sigma_x_k1, Multiply(x_k1, -1.f)); // sigma_x_k1 - x_k1
    std::vector<std::vector<float>> diff_T = Transpose(diff); // Transpose(sigma_x_k1 - x_k1)

    // Kovaryans hesaplama
    for (int i = 0; i < sigma_weights.size(); i++) {
        for (int j = 0; j < p_k1.size(); j++) {
            p_k1[j][0] += sigma_weights[i][1] * (sigma_p[0][i] - x_k1[0][i]) * diff_T[i][0];
            p_k1[j][1] += sigma_weights[i][1] * (sigma_p[1][i] - x_k1[1][i]) * diff_T[i][1];
            p_k1[j][2] += sigma_weights[i][1] * (sigma_p[2][i] - x_k1[2][i]) * diff_T[i][2];
            p_k1[j][3] += sigma_weights[i][1] * (sigma_p[3][i] - x_k1[3][i]) * diff_T[i][3];
        }
    }
    x_k1 = Daraltma(x_k1);

    return std::make_tuple(x_k1, p_k1);
}

//Sgma Bazlý Kalman Filtresi Ýnovasyon
std::tuple<std::vector<std::vector<float>>, std::vector<std::vector<float>>> Sigma_Kalman_Filter_Innovation(std::vector<std::vector<float>> x_k, std::vector<std::vector<float>> p_k, float z_x, float z_vx, float z_y, float z_vy) {
    //Eðer sensörler yalnýzca konum veya hýz ölçüyorsa u(k+1) yer almaz. 
   //Fakat ölçülen deðer doðrudan giriþten (motor gücü ölçümü) etkileniyorsa u(k+1) dahil edilir.
    std::vector<std::vector<float>> z_k = { {z_x},{z_vx},{z_y},{z_vy} }; //Sensörden gelen veriler. (u(k+1) dahil edilmemiþtir.)
    std::vector<std::vector<float>> z_k1(z_k.size(), std::vector<float>(z_k[0].size(), 0)); // 4x1 matris
    std::vector<std::vector<float>> sigma_p = Sigma_Points(x_k, p_k);                       // 4 x 9 matris
    std::vector<std::vector<float>> sigma_weights = Agýrlýk_Hesaplama(sigma_p.size());      // 9 x 2 matris

    //Sensörlerden gelen verileri güncelliyoruz.(k + 1 zamanýndaki)
    for (int i = 0; i < sigma_weights.size(); i++) {
        z_k1[0][0] += sigma_weights[i][0] * sigma_p[0][0];
        z_k1[1][0] += sigma_weights[i][0] * sigma_p[1][0];
        z_k1[2][0] += sigma_weights[i][0] * sigma_p[2][0];
        z_k1[3][0] += sigma_weights[i][0] * sigma_p[3][0];
    }

    //Sensörlerden gelen hata kovaryansýný güncelliyoruz.(k + 1 zamanýndaki)
    //P_yy
    std::vector<std::vector<float>> zp_yy(z_k.size(), std::vector<float>(z_k.size(), 0));                 // 4 x 4 matris
    std::vector<std::vector<float>> p_yy_k1_transpose = Transpose(Matrix_Sum(z_k, Multiply(z_k1, -1.f))); // 1 x 4 matris
    for (int i = 0; i < sigma_weights.size(); i++) {
        for (int j = 0; j < z_k.size(); j++) {
            zp_yy[j][0] += sigma_weights[i][1] * (z_k[0][0] - z_k1[0][0]) * p_yy_k1_transpose[0][0];
            zp_yy[j][1] += sigma_weights[i][1] * (z_k[1][0] - z_k1[1][0]) * p_yy_k1_transpose[0][1];
            zp_yy[j][2] += sigma_weights[i][1] * (z_k[2][0] - z_k1[2][0]) * p_yy_k1_transpose[0][2];
            zp_yy[j][3] += sigma_weights[i][1] * (z_k[3][0] - z_k1[3][0]) * p_yy_k1_transpose[0][3];
        }
    }

    //P_xy
    // 4 x 1 -> 4 x 9
    x_k = Genisletme(x_k);
    std::vector<std::vector<float>> zp_xy(z_k.size(), std::vector<float>(z_k.size(), 0)); // 4 x 4 matris
    std::vector<std::vector<float>> p_xy_k1_transpose = p_yy_k1_transpose;                // 1 x 4 matris
    for (int i = 0; i < sigma_weights.size(); i++) {
        for (int j = 0; j < x_k.size(); j++) {
            zp_xy[j][0] += sigma_weights[i][1] * (sigma_p[0][i] - x_k[0][i]) * p_xy_k1_transpose[0][0];
            zp_xy[j][1] += sigma_weights[i][1] * (sigma_p[1][i] - x_k[1][i]) * p_xy_k1_transpose[0][1];
            zp_xy[j][2] += sigma_weights[i][1] * (sigma_p[2][i] - x_k[2][i]) * p_xy_k1_transpose[0][2];
            zp_xy[j][3] += sigma_weights[i][1] * (sigma_p[3][i] - x_k[3][i]) * p_xy_k1_transpose[0][3];
        }
    }

    //x_k1 4x9 -> 4x1
    x_k = Daraltma(x_k);

    float lambda = 1e-3; // Küçük bir pozitif deðer
    for (int i = 0; i < zp_yy.size(); i++) {
        zp_yy[i][i] += lambda;
    }
    
    //Kalman kazancý
    std::vector<std::vector<float>> kalman_kaz = Matrix_Multiply(zp_xy, Mat_4x4_Ters(zp_yy)); // 4 x 4 matris

    //Matrisleri güncelleþtirme
    x_k = Matrix_Sum(x_k, Matrix_Multiply(kalman_kaz, Matrix_Sum(z_k, Multiply(z_k1, -1.f))));
    p_k = Matrix_Sum(p_k, Multiply(Matrix_Multiply(kalman_kaz, Matrix_Multiply(zp_yy, Transpose(kalman_kaz))), -1.f));

    return std::make_tuple(x_k, p_k);
}

//Klasik Kalman Filtresi için tahmin fonksiyonu
std::tuple<std::vector<std::vector<float>>, std::vector<std::vector<float>>, std::vector<std::vector<float>>, std::vector<std::vector<float>>> Kalman_Filter_Predict(float delta_t, float x_k, float v_x, float y_k, float v_y,
    float a_x, float a_y) {
    float x_k1{}, v_xk1{}, y_k1{}, v_yk1{}, z_k1{}, v_zk1{};
    float sigma_x2 = 1, sigma_xv2 = 1, cov_xv = 1;//rastgele belirlenmiþtir.
    float sigma_y2 = 1, sigma_yv2 = 1, cov_yv = 1;


    float sigma_x2_k1{}, sigma_v2x_k1{}, cov_xv_k1{};//X
    float sigma_y2_k1{}, sigma_v2y_k1{}, cov_yv_k1{};//Y

    std::vector<std::vector<float>> a = { {a_x,a_y,0,0} };//Ývme Vektörü

    std::vector<std::vector<float>> A = { {1,delta_t,0,0},{0,1,0,0},{0,0,1,delta_t},{0,0,0,1} };//Sistem Matrisi

    std::vector<std::vector<float>> B = { {0.5f * delta_t * delta_t},{delta_t},{0.5f * delta_t * delta_t},{delta_t} };//Kontrol Matrisi

    std::vector<std::vector<float>> Q = { {1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1} };//Gürültü Kovaryans Matrisi(Birim matris seçilmiþtir, deðiþtirilebilir.)

    std::vector<std::vector<float>> matrix_k = { {x_k},{v_x},{y_k},{v_y} };//Durum Matrisi(k zamanýndaki)

    std::vector<std::vector<float>> matrix_k1 = { {x_k1},{v_xk1},{y_k1},{v_yk1} };

    std::vector<std::vector<float>> matrix_pk = { {sigma_x2,cov_xv,0,0},{cov_xv,sigma_xv2,0,0},{0,0,sigma_y2,cov_yv},{0,0,cov_yv,sigma_yv2} };//Hata Kovaryans Matrisi(k zamanýndaki)

    std::vector<std::vector<float>> matrix_pk1 = { {sigma_x2_k1,cov_xv_k1,0,0},{cov_xv_k1,sigma_v2x_k1,0,0},{0,0,sigma_y2_k1,cov_yv_k1},{0,0,cov_yv_k1,sigma_v2y_k1} };//Hata Kovaryans Matrisi(k + 1 zamanýndaki)

    //Durum Tahmini
    matrix_k1 = Matrix_Sum(Matrix_Multiply(A, matrix_k), Matrix_Multiply(B, a));

    //Hata Kovaryansý Tahmini
    matrix_pk1 = Matrix_Sum(Matrix_Multiply(A, Matrix_Multiply(matrix_pk, Transpose(A))), Q);

    //Deðerleri bilinmeyenlere atama
    x_k1 = matrix_k1[0][0]; v_xk1 = matrix_k1[1][0]; sigma_x2_k1 = matrix_pk1[0][0]; cov_xv_k1 = matrix_pk1[0][1]; sigma_v2x_k1 = matrix_pk1[1][1]; //X ekseni verilerinin güncellenmesi
    y_k1 = matrix_k1[2][0]; v_yk1 = matrix_k1[3][0]; sigma_y2_k1 = matrix_pk1[2][2]; cov_yv_k1 = matrix_pk1[2][3]; sigma_v2y_k1 = matrix_pk1[3][3]; //Y ekseni verilerinin güncellenmesi

    return std::make_tuple(matrix_k, matrix_k1, matrix_pk, matrix_pk1);
}

//Klasik Kalman Filtresi için inovasyon fonksiyonu
std::tuple<std::vector<std::vector<float>>, std::vector<std::vector<float>>> Kalman_Filter_Innovation(std::vector<std::vector<float>> x_k, std::vector<std::vector<float>> p_k, float z_x, float z_vx, float z_y, float z_vy) {
    float sigma_zx2 = 1, sigma_zvx2 = 1, cov_zxv = 1;//X ekseni için rastgele belirlenmiþtir.
    float sigma_zy2 = 1, sigma_zyv2 = 1, cov_zyv = 1;//Y ekseni için rastgele belirlenmiþtir.

    std::vector<std::vector<float>> H = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };//Gözlem Matrisi

    std::vector<std::vector<float>> R = { {sigma_zx2,cov_zxv,0,0},{cov_zxv,sigma_zvx2,0,0},{0,0,sigma_zy2,cov_zyv},{0,0,cov_zyv,sigma_zyv2,0,0} };//Gözlem Gürültüsü Kovaryans Matrisi

    std::vector<std::vector<float>> z_k = { {z_x},{z_vx},{z_y},{z_vy} };//Gözlem Vektörü

    std::vector<std::vector<float>> brm_mat = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };//Birim matris

    //Kalman Kazanç denklemi
    std::vector<std::vector<float>> Kal_Kaz = Matrix_Multiply(Matrix_Multiply(p_k, Transpose(H)), Mat_4x4_Ters(Matrix_Sum(Matrix_Multiply(H, Matrix_Multiply(p_k, Transpose(H))), R)));

    //Durum tahmin matrisi güncellenir
    x_k = Matrix_Sum(x_k, Matrix_Multiply(Kal_Kaz, Matrix_Sum(z_k, Multiply(Matrix_Multiply(H, x_k), -1.f))));

    //Hata Kovaryans matrisi güncellenir
    p_k = Matrix_Multiply(Matrix_Sum(brm_mat, Multiply(Matrix_Multiply(Kal_Kaz, H), -1.f)), p_k);

    return std::make_tuple(x_k, p_k);
}


int main() {
    sf::RenderWindow window(sf::VideoMode(800, 600), "Kalman Filter");

    //Gerçek
    std::vector<sf::Vector2f> pos;

    //Kalman Tahmini 
	int v_x = 100, v_y = 100;
    int a_x = 100, a_y = 100;
    int t = 0;

    std::vector<std::vector<float>> matrix_k (4,std::vector<float>(1,0));
    std::vector<std::vector<float>> matrix_k1(4, std::vector<float>(1, 0));
    std::vector<std::vector<float>> matrix_pk(4, std::vector<float>(4, 0));
    std::vector<std::vector<float>> matrix_pk1(4, std::vector<float>(4, 0));

    std::vector<sf::Vector2f> Kalman_pos;

    sf::CircleShape cursor;
    cursor.setFillColor(sf::Color::Red);
    
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        // Pencereye göre fare konumu alýnýr.
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);

        pos.push_back(cursor.getPosition());

        cursor.setPosition(static_cast<float>(mousePos.x), static_cast<float>(mousePos.y));

        std::tie(matrix_k1, matrix_pk1) = Sigma_Kalman_Filter_Prediction(t, static_cast<float>(mousePos.x), v_x, static_cast<float>(mousePos.y), v_y, a_x, a_y);
        std::tie(matrix_k, matrix_pk) = Sigma_Kalman_Filter_Innovation(matrix_k1, matrix_pk1, mousePos.x, v_x, mousePos.y, v_y);

        Kalman_pos.push_back(sf::Vector2f(matrix_k[0][0] + 10, matrix_k[2][0] + 10));

        window.clear();

        for (const auto& pos_i : pos) {
            sf::CircleShape circle(3);
            circle.setPosition(pos_i);
            circle.setFillColor(sf::Color::Green);
            window.draw(circle);
        }

        for (const auto& pos_i : Kalman_pos) {
            sf::CircleShape circle(3);
            circle.setPosition(pos_i);
            circle.setFillColor(sf::Color::Blue);
            window.draw(circle);
        }

        t += 1.0 / 60.0;

        std::cout << std::endl << "Zaman: " << t << std::endl;
        std::cout << "Gercek Konum:    " << mousePos.x << "    " << mousePos.y << std::endl;
		std::cout << "Kalman Tahmini: " << matrix_k[0][0] << "     " << matrix_k[2][0] << std::endl;
		std::cout << std::endl;
		std::cout << "X Ekseninde Aralarindaki Fark: " << (mousePos.x - matrix_k[0][0]) << std::endl;
		std::cout << "Y Ekseninde Aralarindaki Fark: " << (mousePos.y - matrix_k[2][0]) << std::endl;

        window.display();
    } 
   
    return 0;
}