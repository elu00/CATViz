#include <fstream>
#include <string>
#include <sstream>
#include "glm/vec2.hpp"
#include "glm/geometric.hpp"

using glm::vec2;
using std::string;
using std::endl;

class CAT
{
    public:
        vec2 i, j, k;
        double a_ij, a_jk, a_ki;
        string to_string()
        {
            // header
            std::stringstream ss;
            ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
                << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
                << "<svg width=\"500\" height=\"500\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
            /*
               ss  << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl
               << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"" << endl
               << " \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl 
               << "<svg width=\"" << 100 
               << "\" height=\"" << 100 << "\">" << endl
               << "<title>Collection of Geometry</title>" << endl;
               */
            // points
            ss << "<circle cx=\"" << i.x << "\" cy=\"" << i.y << "\" r=\"5\"/>" << endl;
            ss << "<circle cx=\"" << j.x << "\" cy=\"" << j.y << "\" r=\"5\"/>" << endl;
            ss << "<circle cx=\"" << k.x << "\" cy=\"" << k.y << "\" r=\"5\"/>" << endl;
            // arcs
            ss << to_arc(i, j, a_ij) << endl << to_arc(j, k, a_jk) << endl << to_arc(k, i, a_ki) << endl;

            // footer
            ss << "</svg>";

            return ss.str();

        }
        //overload for just sidelengths
        CAT (double ij, double jk, double ki, double a1, double a2, double a3): a_ij(a1), a_jk(a2), a_ki(a3)
    {
    }
        // overload for specifying coordinates
        CAT (double i1, double i2, double j1, double j2, double k1, double k2, double a1, double a2, double a3): a_ij(a1), a_jk(a2), a_ki(a3)
    {
        i = vec2(i1, i2);
        j = vec2(j1, j2);
        k = vec2(k1, k2);
    }
        void to_svg(string filename)
        {
            std::ofstream f (filename, std::ofstream::out);
            f << to_string();
            f.close();
            return;
        }

    private:
        string to_arc(vec2 a, vec2 b, double angle)
        {
            double radius = glm::distance(a, b) / (2 * angle);
            string largeArcFlag = std::abs(2*angle) <= 3.14159265358979323846264 ? "0" : "1";
            // sweep flag is 0 if going outward, 1 if going inward
            string sweepFlag = angle < 0 ? "1" : "0";
            std::stringstream ss;
            ss << "<path d=\"M" << a.x << "," << a.y << " A" << radius << "," 
                << radius << " 0 " <<  largeArcFlag << " " << sweepFlag << " " 
                << b.x << "," << b.y << "\" fill=\"red\" stroke=\"blue\" stroke-width=\"5\" />"; 

            return ss.str();

        }
        vec2 center(vec2 a, vec2 b, double angle)
        {
            // |b-a| = 2 alpha r, so radius = |b-a| / 2aalpha
            double radius = glm::distance(a, b) / (2 * angle);
            vec2 diff = b - a;
            // counterclockwise rotation by 90 degrees
            vec2 offset = vec2(-diff.y, diff.x);
            if (angle <= 0)
            {
                offset *= radius/glm::length(diff);
            }
            else
            {
                offset *= -radius/glm::length(diff);
            }
            return (a + b)/2.f + offset;
        }

};


