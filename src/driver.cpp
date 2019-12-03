#include <iostream>
#include "util.cpp"

using std::cout;
using std::cin;

int main()
{
    /*
    double i, j, k, a, b, c;
    string filename;
    cout << "Enter i" << endl;
    cin >> i;
    cout << "Enter j" << endl;
    cin >> j;
    cout << "Enter k" << endl;
    cin >> k;
    CAT triangle = CAT(i, j, k, a, b, c, true);
    cout << i << endl << j << endl << k << endl << a << endl << b << endl << c << endl;
    */
    double i1, i2, j1, j2, k1, k2, a, b, c;
    string filename;
    cout << "Enter i1" << endl;
    cin >> i1;
    cout << "Enter i2" << endl;
    cin >> i2;
    cout << "Enter j1" << endl;
    cin >> j1;
    cout << "Enter j2" << endl;
    cin >> j2;
    cout << "Enter k1" << endl;
    cin >> k1;
    cout << "Enter k2" << endl;
    cin >> k2;
    cout << "Coordinates are: " 
        << "(" << i1 << "," << i2 << ")"
        << "(" << j1 << "," << j2 << ")"
        << "(" << k1 << "," << k2 << ")" << endl;

    cout << "Enter alpha" << endl;
    cin >> a;
    cout << "Enter beta" << endl;
    cin >> b;
    cout << "Enter gamma" << endl;
    cin >> c;
    cout << "Enter filename" << endl;
    cin >> filename;
    CAT triangle = CAT(i1, i2, j1, j2, k1, k2, a, b, c, true);
    triangle.to_svg(filename);
    return 1;
}
