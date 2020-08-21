#include <iostream>
#include "util.cpp"
#include "glm/gtx/string_cast.hpp"

using std::cout;
using std::cin;
void SVGTests()
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
    return;
}
void BaryTests() {
    //assert(abs(l2DistSquared(5, 5 * sqrt(2), 5, 0, 0, 0, 1, 0, 0, 0, 1, 0) - 25.0) < 1e-8);
    //assert(abs(l2DistSquared(5, 5 * sqrt(2), 5, 0, 0, 0, 0, 1, 0, 0.5, 0.5, 0) - 6.25) < 1e-8);
    //cout << sqrt(l2DistSquared(5, 5 * sqrt(2), 5, -1e-2, -1e-2, -1e-2, 0.3, 0.7, 0, 0.3, 0.2, 0.5)) << endl;
    //cout << l2DistSquared(5 * sqrt(2), 5, 5, -1e-2, -1e-2, -1e-2, 0.7, 0, 0.3, 0.2, 0.5, 0.3) << endl;
    cout << l2DistSquared(12.322, 2.6971, 10.4546, 0.417028, 0.423979, 0.26072, 0, 0, 0, 0, 0, 0) << endl;
}
int main() {
    //SVGTests();
    CAT triangle = CAT(200, 140, 120, 0.1, -0.6, -0.2, true);
    triangle.to_svg("orig.svg");
    std::ofstream f ("tri.svg", std::ofstream::out);
    f << triangle.triangulation_svg(5);
    f.close();
    vector<double> finalPoints;
    vector<int> triangles;
    std::tie(finalPoints, triangles) = triangle.triangulation(5);
    for (int i = 0; i < finalPoints.size()/2; i++) {
        cout << glm::to_string(triangle.planeToBary(vec2(finalPoints[2*i], finalPoints[2*i+1]))) << endl;
    }
    return 0;
}


