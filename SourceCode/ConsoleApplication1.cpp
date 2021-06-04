#include "function.cpp"
#include "fstream"
#include <string>
int main()
{
    int i = 0;
    int j;
    kvz_pixel my_ref_left [17] = { 160, 150, 140, 130, 120, 110, 50, 40, 30, 20, 20, 20, 20, 20, 20, 20 };
    kvz_pixel my_ref_top[17] = { 160, 150, 140, 130, 50, 40, 40,  30, 20, 20, 20, 20, 20, 20, 20, 20 };
    kvz_intra_references refs_u;
    for (i = 0; i <= 16; i++) {
        refs_u.ref.left[i] = my_ref_left[i];
        refs_u.ref.top[i] = my_ref_top[i];
        refs_u.filtered_ref.left[i] = my_ref_left[i];
        refs_u.filtered_ref.top[i] = my_ref_top[i];
        refs_u.filtered_initialized = true;
    }
    
    for (j = 0; j <= 34 ; j++) {

        kvz_pixel* pred = new kvz_pixel[32 * 32 + 32];
        kvz_intra_predict(&refs_u, 3 /*log2_width*/, j /*mode*/, COLOR_U, pred, true);
        ofstream file;
        file.open(std::to_string(j) + ".txt");
        for (i = 0; i <= 63; i++) {
            cout << unsigned(pred[i]) << "\n";
            file << unsigned(pred[i]) << "\n";
        }
        file.close();
    }
    
}

