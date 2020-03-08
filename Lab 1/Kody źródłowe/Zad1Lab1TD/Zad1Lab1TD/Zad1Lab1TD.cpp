#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    float a = 9;
    float b = 8;
    float c = 2;

    float det = b * b - 4 * a * c;
    cout << "Wyznacznik: " << det << "." << endl;
    if (det < 0)
    {
        cout << "Brak miejsc zerowych." << endl;
    }
    else if(det == 0)
    {
        float t = -b / 2 * a;

        cout << "Istnieje jedno miejsce zerowe t=" << t << "." << endl;
    }
    else
    {
        float t1 = (-b + sqrt(det)) / 2 * a;
        float t2 = (-b - sqrt(det)) / 2 * a;

        cout << "Istnieją dwa miejsca zerowe t1=" << t1 << ", t2=" << t2 << "." << endl;
    }

    //Do wykresu*******************
    ofstream save("data.txt");
    
    for (double i = -10; i < 10; i = i + 0.01)
    {
        //cout << i << endl;
        float fun = a * i * i + b * i + c;
        //cout << fun << endl;
        save << fun << endl;
    }

    save.close();
}