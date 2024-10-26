#include <bits/stdc++.h>
typedef long long ll;
#define forr(i, j, k) for (ll i = j; i < k; i++)
using namespace std;
void gauss_siedel(){
    int n;
    cout << "Enter the number of variables/equations: ";
    cin >> n;

    // Vectors to hold coefficients and constants
    vector<vector<float>> a(n, vector<float>(n)); // Coefficients
    vector<float> r(n);                           // Constants
    vector<float> x(n, 0);                        // Solution vector (initial guesses)
    vector<float> x_new(n);                       // Updated solution
    float e;                                      // Precision

    
    cout << "Enter the coefficients for each equation (row-wise):\n";
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cin >> a[i][j]; 
        }
        cin >> r[i]; 
    }

    cout << setprecision(6) << fixed;
    cout << "Enter precision: ";
    cin >> e;

    int step = 1;
    float max_error; 
    do
    {
        
        for (int i = 0; i < n; i++)
        {
            float sum = r[i];
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum -= a[i][j] * x[j];
                }
            }
            x_new[i] = sum / a[i][i]; 
        }

        
        cout << step << "\t";
        for (int i = 0; i < n; i++)
        {
            cout << "x" << i + 1 << "= " << x_new[i] << ", ";
        }
        cout << endl;

        
        max_error = 0; 
        for (int i = 0; i < n; i++)
        {
            max_error = max(max_error, fabs(x_new[i] - x[i])); 
        }

        
        x = x_new; 

        step++;
    } while (max_error > e); 
}
int main()
{
    gauss_siedel();
    return 0;
}
