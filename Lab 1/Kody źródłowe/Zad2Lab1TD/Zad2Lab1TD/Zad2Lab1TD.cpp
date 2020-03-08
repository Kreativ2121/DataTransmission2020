#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main()
{
    int a = 9;
    int b = 8;
    int c = 2;

    float pi = 3.14159265359;

    //Funkcje y, z, u

    ofstream savey("data_y.txt");
    ofstream savez("data_z.txt");
    ofstream saveu("data_u.txt");
    ofstream saveOX("data_OX.txt");
    for (double i = 0; i <= 1; i = i + 1./22050)
    {
        //cout << i << endl;
        float x = a * i * i + b * i + c;
        float y = 2 * x * x + 12 * cos(i);
        savey << y << endl;
        float z = sin(2 * pi * 7 * i) * x - 0.2 * log10(abs(y) + pi);
        savez << z << endl;
        float u = sqrt(abs(y * y * z)) - 1.8 * sin(0.4 * i * z * x);
        saveu << u << endl;
        saveOX << i << endl;
    }

    savey.close();
    savez.close();
    saveu.close();
    saveOX.close();

    //Funkcja v

    ofstream savev("data_v.txt");

    for (double i = 0; i <= 1; i = i + 1./22050)
    {
        float v;
        if (i < 0.22)
            v = (1 - 7 * i) * sin((2 * pi * i * 10) / (i + 0.04));
        else if (i < 0.7)
            v = 0.63 * i * sin(125 * i);
        else
            v = pow(i, -0.662) + 0.77 * sin(8 * i);

        savev << v << endl;
    }

    savev.close();

    //Funkcja p

    ofstream savep("data_p.txt");

    //int N = 2;
    //int N = 4;
    int N = 98;

    for (double i = 0; i <= 1; i = i + 1. / 22050)
    {
        float p = 0;
        for (int n = 1; n < N; n++)
        {
            p +=(cos(12 * i * n * n) + cos(16 * i * n)) / (n * n);
        }
        savep << p << endl;
    }

    savep.close();
}
