#include <iostream>
#include "svg_gen.hpp"

using std::cout;
using std::cin;

int main()
{
    cout << "Enter i1" << endl;
    double i1;
    while (!(cin >> i1))
    {
        cout << "Please enter a valid input";
    }
    cout << "Enter i2" << endl;
    double i2;
    while (!(cin >> i2))
    {
        cout << "Please enter a valid input";
    }
    cout << "Enter j1" << endl;
    double j1;
    while (!(cin >> j1))
    {
        cout << "Please enter a valid input";
    }
    cout << "Enter j2" << endl;
    double j2;
    while (!(cin >> j2))
    {
        cout << "Please enter a valid input";
    }
    cout << "Enter k1" << endl;
    double k1;
    while (!(cin >> k1))
    {
        cout << "Please enter a valid input";
    }
    cout << "Enter k2" << endl;
    double k2;
    while (!(cin >> k2))
    {
        cout << "Please enter a valid input";
    }
    cout << "Coordinates are: " 
            << "(" << i1 << "," << i2 << ")"
            << "(" << j1 << "," << j2 << ")"
            << "(" << k1 << "," << k2 << ")" << endl;
    cout << "Enter filename" << endl;
    string filename;




}
