#include <fstream>
#include <string>
#include <sstream>
#include "glm/vec2.hpp"

using glm::vec2;
using std::string;
using std::endl;

class CAT
{
    public:
        vec2 i, j, k;
        double a_ij, a_jk, a_ki;
        CAT (vec2 i, vec2 j, vec2 k, double a1, double a2, double a3);
        string to_string();
        void to_svg(string filename);
    private:
        vec2 center(vec2 a, vec2 b, double angle);


};
CAT::CAT (vec2 i, vec2 j, vec2 k, double a1, double a2, double a3): a_ij(a1), a_jk(a2), a_ki(a3)
{
    // 
}
string CAT::to_string()
{
    std::stringstream ss;
    ss  << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"" << endl
        << " \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl 
        << "<svg width=\"" << 100 
        << "\" height=\"" << 100 << "\">" << endl
        << "<title>Bar</title>" << endl;


    return ss.str();
}

void CAT::to_svg(string filename)
{
    std::ofstream f (filename, std::ofstream::out);
    f << to_string();
    f.close();
    return;
}


