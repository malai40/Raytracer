// Raytracer (reads NFF file format and exports PNG images)
// Source on NFF file format: https://paulbourke.net/dataformats/nff/nff1.html

#include <istream>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "LodePNG/lodepng.h"
#include "Eigen/Dense"

// Declare global r, g, b
float r = 1;
float g = 1;
float b = 1;

// Declare shading parameters
// Diffuse component
float kd = 1;
// Specular component
float ks = 0;
// Phong cosine power for highlights
float shine = 0;
// Transmittance (fraction of contribution of the transmitting ray)
float t = 0;
// Index of refraction
float index_ref = 1;

// Declare global background r, g, b
Eigen::Vector4d bg_c(0, 0, 0, 1);

// Declare width and height of final image
unsigned width = 0;
unsigned height = 0;

// Declare the distance from at to eye
float dist_at_to_eye = 0;

// Declare the eye
// Default is (0,0,0), which can be changed with "eye" call
Eigen::Vector3d eye(0, 0, 0);

// Declare the forward
// Default is (0,0,-1), which can be changed with "forward" call
Eigen::Vector3d forward(0, 0, -1);

// Declare the up
// Default is (0,1,0), which can be changed with "up" call
Eigen::Vector3d up(0, 1, 0);

// Declare the angle/zoom
// Stored in radians (radians for calculations)
double fov_rad = 0;

// Declare hither (near plane)
float hither = 0;

// Global vertex array for triangles
std::vector<Eigen::Vector3d> vertices;

// For normals
bool read_normals = false;
// Have a current normal at all times
// Every time a xyz is read in, add the current normal to normals
Eigen::Vector3d normal_global;

// Global normal array for triangles
std::vector<Eigen::Vector3d> normals;

// Read in text instructions to create image file
// Ignores lines starting with the "#" character and blank lines
std::vector<std::string> text_file_to_lines(const char* filename) {
    // Declare the array of lines
    std::vector<std::string> lines;

    // Start reading the file
    std::ifstream f;
    f.open(filename);

    // Use getline to read each line
    std::string line;
    while (std::getline(f, line)) {
        if (line.size() > 0) { //Line is not just return character(s)
            if (line[0] != '#') {
                lines.push_back(line);
            }
        }
    }
    
    // Close the file stream
    f.close();

    // Return lines
    return lines;

}

// Take a string with too many spaces, and return it with just single spaces
std::string remove_extra_whitespace(std::string str) {
    std::string str_final;

    unsigned start = 0;
    unsigned end = 0;
    unsigned i = 0;
    unsigned space_counter = 0;
    //unsigned len = str.length();
    unsigned chars = 0;
  
    //while (i < len)
    //for (int i = 0; i < len; i++) {
    for (int i = 0; i < str.length(); i++) {
        if (str[i] != ' ') {
            chars++;
            
            //if (i == len-1) {
            if (i == str.length()-1) {
                if (chars > 0) {
                    //end = len;
                    end = str.length();
                    //str_final += str.substr(start, end);
                    str_final = str_final + str.substr(start, end);
                }
                //i++; // TODO This is useless since we're at end of string, right?
            }     
            space_counter = 0;
            //i++;
        } else {
            space_counter++;
            if (space_counter == 1) {                
                end = i + 1 - start;
                //str_final += str.substr(start, end);
                str_final = str_final + str.substr(start, end);
                start = i + 1;
                chars = 0;
            } else if (space_counter > 1) {
                start++;
            }
            //i++;
        } 
    }

    // Remove space at beginning and end
    if (str_final[0] == ' ') {
        str_final = str_final.substr(1,str_final.length()-1);
    }
    if (str_final[str_final.length()-1] == ' ') {
        str_final = str_final.substr(0,str_final.length()-1);
    }

    return str_final;
}

// Transform a line of tokens into an array of tokens
std::vector<std::string> line_to_tokens(std::string line0) {
    // Break up a string into tokens by spaces
    std::string line = remove_extra_whitespace(line0);
    
    std::vector<std::string> tokens;

    std::istringstream line_stream (line);

    std::string token;
    while(std::getline(line_stream, token, ' ')) {
        tokens.push_back(token);
    }

    return tokens;

}

// Change the camera's eye location
void get_eye_info(std::vector<std::string> line_tokens) {
    float x = stof(line_tokens[1]);
    float y = stof(line_tokens[2]);
    float z = stof(line_tokens[3]);

    eye(0) = x;
    eye(1) = y;
    eye(2) = z;
}

// Calculate forward vector using "eye"/"from" and "at"
void get_at_info(std::vector<std::string> line_tokens) {
    float x = stof(line_tokens[1]);
    float y = stof(line_tokens[2]);
    float z = stof(line_tokens[3]);
    
    // Calculate vector from eye to at
    float vec_x = x - eye(0);
    float vec_y = y - eye(1);
    float vec_z = z - eye(2);
    
    // Make new vector
    Eigen::Vector3d forward_pre(vec_x, vec_y, vec_z);
    // Save length of this vector as distance from at to eye
    dist_at_to_eye = forward_pre.norm();
    // Normalize
    forward_pre.normalize();
    
    // Set forward to this newly normalized vector
    forward(0) = forward_pre(0);
    forward(1) = forward_pre(1);
    forward(2) = forward_pre(2);   
    
}

// Change the camera's up vector
void get_up_info(std::vector<std::string> line_tokens) {
    float x = stof(line_tokens[1]);
    float y = stof(line_tokens[2]);
    float z = stof(line_tokens[3]);

    up(0) = x;
    up(1) = y;
    up(2) = z;
}

// Change the background color of final image
void get_bg_info(std::vector<std::string> line_tokens) {
    float bg_r_temp = stof(line_tokens[1]);
    float bg_g_temp = stof(line_tokens[2]);
    float bg_b_temp = stof(line_tokens[3]);

    bg_c(0) = bg_r_temp;
    bg_c(1) = bg_g_temp;
    bg_c(2) = bg_b_temp;
    bg_c(3) = 1;
}

// Change the camera's field of view (FOV) using "angle"
// Implement the FOV into the forward
void get_angle_info(std::vector<std::string> line_tokens) {
    float angle_degrees = stof(line_tokens[1]);
    
    // Store as radians in fov
    fov_rad = angle_degrees * (M_PI/180);


}

// Get near plane (hither)
void get_hither_info(std::vector<std::string> line_tokens) {
    float near = stof(line_tokens[1]);
    
    hither = near;

}

// Get color and shading parameters
void color(float &r, float &g, float &b, std::vector<std::string> line_tokens) {
    // r, g, b are floating-point values
    r = stof(line_tokens[1]);
    g = stof(line_tokens[2]);
    b = stof(line_tokens[3]);
    
    // Take in other values
    kd = stof(line_tokens[4]);
    ks = stof(line_tokens[5]);
    shine = stof(line_tokens[6]);
    t = stof(line_tokens[7]);
    index_ref = stof(line_tokens[8]);
 
}

// Get image resolution
void get_resolution_info(std::vector<std::string> line_tokens) {
    float w = stof(line_tokens[1]);
    float h = stof(line_tokens[2]);
    
    width = w;
    height = h;

}

// Add a directional light (sun) to the list of objects to be rendered
std::vector<double> load_sun(std::vector<std::string> line_tokens) {
    float x_td = stof(line_tokens[1]);
    float y_td = stof(line_tokens[2]);
    float z_td = stof(line_tokens[3]);

    std::vector<double> scene_obj;
    // Store a number to indicate shape type
    // 0 is sphere, 1 is sun, 2 is plane, others... TODO
    scene_obj.push_back(1);
    // Then the coordinates of the sphere
    
    scene_obj.push_back(x_td); //x
    scene_obj.push_back(y_td); //y
    scene_obj.push_back(z_td); //z
    // Then the RGB (TODO and other things like shininess)
    // If RGB was supplied, set that too TODO
    if (line_tokens.size() > 4) {
        float r_td = stof(line_tokens[4]);
        float g_td = stof(line_tokens[5]);
        float b_td = stof(line_tokens[6]);
        scene_obj.push_back(r_td);
        scene_obj.push_back(g_td);
        scene_obj.push_back(b_td);
    } else {
        scene_obj.push_back(1); // White light
        scene_obj.push_back(1); // White light
        scene_obj.push_back(1); // White light
    }
   
    scene_obj.push_back(1); //Current alpha //TODO How to implement alpha for optional things?

    return scene_obj;
}

// Add a positional light (bulb) to the list of objects to be rendered
std::vector<double> load_bulb(std::vector<std::string> line_tokens) {
    float x_td = stof(line_tokens[1]);
    float y_td = stof(line_tokens[2]);
    float z_td = stof(line_tokens[3]);

    std::vector<double> scene_obj;
    // Store a number to indicate shape type
    // 0 is sphere, 1 is sun, 2 is plane, others... TODO
    // 4 is bulb (positional light)
    scene_obj.push_back(4);
    // Then the coordinates of the sphere
    
    scene_obj.push_back(x_td); //x
    scene_obj.push_back(y_td); //y
    scene_obj.push_back(z_td); //z
    // Then the RGB (TODO and other things like shininess)
    // If RGB was supplied, set that too TODO
    if (line_tokens.size() > 4) {
        float r_td = stof(line_tokens[4]);
        float g_td = stof(line_tokens[5]);
        float b_td = stof(line_tokens[6]);
        scene_obj.push_back(r_td);
        scene_obj.push_back(g_td);
        scene_obj.push_back(b_td);
    } else {
        scene_obj.push_back(1); // White light
        scene_obj.push_back(1); // White light
        scene_obj.push_back(1); // White light
    }
   
    scene_obj.push_back(1); //Current alpha //TODO How to implement alpha for optional things?

    return scene_obj;
}

// Add a sphere to the list of objects to be rendered
std::vector<double> load_sphere(std::vector<std::string> line_tokens) {
    float x_td = stof(line_tokens[1]);
    float y_td = stof(line_tokens[2]);
    float z_td = stof(line_tokens[3]);
    float radius = stof(line_tokens[4]);
    
    // Create the equation based on the center and radius of the sphere
    // Equation of a sphere:
    // (x-a)^2 + (y-b)^2 + (z-c)^2 = r^2, where (a,b,c) is the sphere center
    // Equation for sphere from https://www.expii.com/t/equation-of-a-sphere-1321

    std::vector<double> scene_obj;
    // Store a number to indicate shape type
    // 0 is sphere, 1 is sun, 2 is plane, others... TODO
    scene_obj.push_back(0);
    // Then the coordinates of the sphere
    
    scene_obj.push_back(x_td); //x
    scene_obj.push_back(y_td); //y
    scene_obj.push_back(z_td); //z
    scene_obj.push_back(radius); //radius
    // Then the RGB (TODO and other things like shininess)
    scene_obj.push_back(r); //Current R
    scene_obj.push_back(g); //Current G
    scene_obj.push_back(b); //Current B
    scene_obj.push_back(1); //Current alpha //TODO How to implement alpha for optional things?

    return scene_obj;
}


// Clamp number to 0
float clamp_float_0(float n) {
    if (n < 0.0) {
        return 0;
    }

    return n;
}

// Clamp number to 1
float clamp_float_1(float n) {
    if (n > 1.0) {
        return 1;
    }

    return n;
}

// Clamp number to [0,1] range
float clamp_float_0_1(float n) {
    n = clamp_float_0(n);
    n = clamp_float_1(n);

    return n;
}

// Convert from linear color to sRGB (gamma) color
float linear_to_srgb(float c) {
    if (c <= 0.0031308) {
        return (12.92 * c);
    } else {
        return (1.055 * std::pow(c,(1/2.4)) - 0.055);
    }
}



unsigned convert_linear_color_to_unsigned_srgb(float val) {
    // r, g, b are floating-point values

    // Clamp <0 to 0 and >1 to 1
    val = clamp_float_0_1(val);

    // Apply a gamma function (convert from linear to sRGB)
    val = linear_to_srgb(val);

    // Scale up linearly to the 0-255 byte range
    unsigned val_f = val * 255;

    return val_f;
}

unsigned convert_linear_color_to_unsigned_linear(float val) {
    // r, g, b are floating-point values

    // Clamp <0 to 0 and >1 to 1
    val = clamp_float_0_1(val);

    // Scale up linearly to the 0-255 byte range
    unsigned val_f = val * 255;

    return val_f;
}


// From LodePNG CPP example (LodePNG's manual)
// https://raw.githubusercontent.com/lvandeve/lodepng/master/examples/example_encode.cpp
//std::vector<unsigned char> color_pixel(std::vector<unsigned char> &image, unsigned width, unsigned x, unsigned y, unsigned red, unsigned green, unsigned blue, unsigned alpha) {
void color_pixel(std::vector<unsigned char> &image, unsigned width, unsigned x, unsigned y, unsigned red, unsigned green, unsigned blue, unsigned alpha) {
    image[((y * width) + x)*4 + 0] = red; //red;
    image[((y * width) + x)*4 + 1] = green; //green;
    image[((y * width) + x)*4 + 2] = blue; //blue;
    image[((y * width) + x)*4 + 3] = alpha; //alpha; // Ranges from 0 to 255, not 0.0 to 1.0

    //return image;
}

// Return all intersections, not just the minimum
// Another function will return the minimum intersection
// TODO Did I use enough of this to bother citing?
// Explanation of how to interpret signs of t's from https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection.html
std::vector<std::vector<double>> get_intersect_all(Eigen::Vector3d origin_v, Eigen::Vector3d direction_v, std::vector<std::vector<double>> scene_objs) {
    double t_min = NAN; //Holds info on where collision happened
    int i_min = -1; // Holds info on which object the minimum collision happened with. If stays -1, no collision with this pixel's ray.
    double t0 = NAN;
    double t1 = NAN;

    std::vector<double> a_solution(2);
    std::vector<std::vector<double>> all_solutions;

    // Break the origin and direction of this pixel down
    // Origin
    double o_x = origin_v(0);
    double o_y = origin_v(1);
    double o_z = origin_v(2);

    // Direction
    double d_x = direction_v(0);
    double d_y = direction_v(1);
    double d_z = direction_v(2);

    // Loop through all of the objects and see if there's an intersection
    // // with this pixel's ray.
    // We need to solve various intersection equations depending on the
    // // object type
    for (int i = 0; i < scene_objs.size(); i++) {
        if (scene_objs[i][0] == 0) {
            // This is a sphere
            // Get the center and radius of this sphere
            double center_x = scene_objs[i][1];
            double center_y = scene_objs[i][2];
            double center_z = scene_objs[i][3];
            double radius = scene_objs[i][4];


            // TODO Translate to vector math to find t.
            double a = pow(d_x, 2) + pow(d_y, 2) + pow(d_z, 2);
            double b = (2 * o_x * d_x) + (2 * o_y * d_y) + (2 * o_z * d_z) - (2 * d_x * center_x) - (2 * d_y * center_y) - (2 * d_z * center_z);
            double c = (pow(o_x, 2) - (2 * o_x * center_x) + pow(center_x, 2) 
                        + pow(o_y, 2) - (2 * o_y * center_y) + pow(center_y, 2) 
                        + pow(o_z, 2) - (2 * o_z * center_z) + pow(center_z, 2) 
                        - pow(radius, 2));
            
            // Use Quadratic formula to solve for t
            // Find t when the operation is subtract
            // Solve for just the sqrt part first. That will tell us if there are any real solutions.
            double dis = sqrt(pow(b, 2) - (4 * a * c));
            // Do nothing if dis is negative
            // If dis == 0, there's just one t.
            // If dis > 0, there are two t's.
            t0 = NAN;
            t1 = NAN;

            

            if (dis < 0) {
                // There are no real solutions
                t0 = NAN;
                t1 = NAN;
            } else if (dis == 0) {
                //printf("Obj #%i, dis is %f\n", i, dis);
                // There is 1 real solution
                t0 = (((-1 * b) - sqrt(pow(b, 2) - (4 * a * c))) / (2 * a));
                a_solution[0] = i;
                a_solution[1] = t0;

                all_solutions.push_back(a_solution);
                
            } else if (dis > 0) {
                //printf("Obj #%i, dis is %f\n", i, dis);
                // There are two real solutions
                // Find t when the operation is add
                t0 = (((-1 * b) - sqrt(pow(b, 2) - (4 * a * c))) / (2 * a));
                // Find t when the operation is add
                t1 = (((-1 * b) + sqrt(pow(b, 2) - (4 * a * c))) / (2 * a));

                //printf("t0 is %f, t1 is %f\n", t0, t1);
                
                a_solution[0] = i;
                a_solution[1] = t0;

                all_solutions.push_back(a_solution);

                //a_solution[0] = i;
                a_solution[1] = t1;

                all_solutions.push_back(a_solution);
            }

        } else if (scene_objs[i][0] == 2 || scene_objs[i][0] == 3) {
            // This is a plane or triangle
            
            // Get the A, B, C, D of this plane
            double a = scene_objs[i][1];
            double b = scene_objs[i][2];
            double c = scene_objs[i][3];
            double d = scene_objs[i][4];

            // Solve this equation: et = -1 * f
            double e = (a * d_x) + (b * d_y) + (c * d_z);
            //double f = -1 * ((a * o_x) + (b * o_y) + (c * o_z) + d);
            double f = 0 - ((a * o_x) + (b * o_y) + (c * o_z) + d);

            double t0 = f / e;

            // If this is a triangle, d_x, d_y, and d_z must be 
            // // within the triangle's vertices.
            // // Set t to NAN if (d_x, d_y, d_z) is outside of the triangle
            if (scene_objs[i][0] == 3 && !std::isnan(t0)) {
                //std::cout << "Found a triangle" << std::endl;
                // This is a triangle
                // Get the vertices of the triangle
                Eigen::Vector3d vert0 = vertices[(int)scene_objs[i][9]];
                Eigen::Vector3d vert1 = vertices[(int)scene_objs[i][10]];
                Eigen::Vector3d vert2 = vertices[(int)scene_objs[i][11]];

                // Turn this point into a "vector"
                Eigen::Vector3d this_point(o_x + d_x*t0, o_y + d_y*t0, o_z + d_z*t0);

                // Get vectors perpendicular to the edges
                Eigen::Vector3d ab = vert0 - vert1;
                Eigen::Vector3d ac = vert0 - vert2;
                Eigen::Vector3d bc = vert1 - vert2;

                // Find the normal
                // Either of the below works
                //Eigen::Vector3d norm_bary = ab.cross(ac);
                Eigen::Vector3d norm_bary(a, b, c);

                Eigen::Vector3d a1 = (vert2 - vert0).cross(norm_bary);
                Eigen::Vector3d a2 = (vert1 - vert0).cross(norm_bary);

                Eigen::Vector3d e1 = a1 / a1.dot(vert1 - vert0);
                Eigen::Vector3d e2 = a2 / a2.dot(vert2 - vert0);

                double w2 = e2.dot(this_point - vert0);
                double w1 = e1.dot(this_point - vert2);
                double w0 = 1 - w1 - w2;
               

/*
                // Get vector perpendicular to the normal and the triangle's edge
                Eigen::Vector3d perp_ab = ab.cross(norm_bary);
                Eigen::Vector3d perp_ac = ac.cross(norm_bary);
                Eigen::Vector3d perp_bc = bc.cross(norm_bary);

                // Dot product the perp by one of the touching edges
                // // and divide by the result.
                perp_ab = perp_ab / (perp_ab.dot(ac));
                perp_ac = perp_ac / (perp_ac.dot(bc));

                double b_ab = perp_ab.dot(this_point - vert0);
                double b_ac = perp_ac.dot(this_point - vert2);
                double b_bc = 1 - b_ab - b_ac;
*/
                // Test: Is perp_ab dot ac == 1?
                /*
                printf("Testing: e2 dot (p2 - p0) is %f\n", e2.dot(vert2 - vert0));
                printf("Testing: e1 dot (p1 - p0) is %f\n", e1.dot(vert1 - vert0));
                printf("Testing: e2 dot (p2 - p1) is %f\n", e2.dot(vert2 - vert1));
                */
                // Get w's by (vertex - this_point), then dot perp
                //double w0 = perp_ab.dot(this_point - vert0);
                //double w1 = perp_ac.dot(this_point - vert0);
                //double w2 = 1 - w0 - w1;

                //double w0 = b0;
                //double w1 = b1;
                //double w2 = b2;
                //double w2 = this_point.dot(perp_bc);

                // If all w's are positive, point is inside the triangle
                if (w0 < 0 || w1 < 0 || w2 < 0) {
                    t0 = NAN;
                /*
                } else {
                    //printf("Point (%f, %f, %f) is within triangle (%f, %f, %f), (%f, %f, %f), (%f, %f, %f)\n", this_point[0], this_point[1], this_point[2], vert0[0], vert0[1], vert0[2], vert1[0], vert1[1], vert1[2], vert2[0], vert2[1], vert2[2]);
                */
                }
                
                // Use barycentric coordinates to see if (d_x, d_y, d_z)
                // // is in this triangle
                // Find w0, w1, and w2
                // If any of the w's are outside of [0,1], set t0 to NAN


            }

            if (!std::isnan(t0)) {
                //std::cout << "Found a triangle to light" << std::endl;
                a_solution[0] = i;
                a_solution[1] = t0;

                all_solutions.push_back(a_solution); 
            }

        }
    }

    return all_solutions;

}

// For shadow calculation, move origin point a millionth of a unit away from true origin before computing
// // to avoid self-shadowing due to touching the surface of the object
Eigen::Vector3d move_origin_for_shadow(Eigen::Vector3d origin, Eigen::Vector3d direction) {
    Eigen::Vector3d new_origin;

    double t = 0.000001;
    new_origin = origin + (t * direction);

    return new_origin;
}

// Get the minimum solution in a list of all of them
std::vector<double> get_min_solution(std::vector<std::vector<double>> solutions) {
    double min_i = -1;
    double min_t = 9999999999; // TODO Lazy.
    
    for (int i = 0; i < solutions.size(); i++) {
        if (min_t > solutions[i][1] && solutions[i][1] > 0) {
            min_t = solutions[i][1];
            min_i = solutions[i][0];
        }
    }

    std::vector<double> final_sol{min_i, min_t};

    return final_sol;
}

// Get the minimum solution for solving shadow rays
// The minimum must not be of the same object that it came from.
std::vector<double> get_min_solution_shadow(std::vector<std::vector<double>> solutions, int this_obj) {
    double min_i = -1;
    double min_t = 9999999999; // TODO Lazy.
    
    for (int i = 0; i < solutions.size(); i++) {
        if (min_t > solutions[i][1] && this_obj != (int)solutions[i][0] && solutions[i][1] >= 0) {
            min_t = solutions[i][1];
            min_i = solutions[i][0];
        }
    }

    std::vector<double> final_sol{min_i, min_t};

    return final_sol;
}

// Add a plane to the list of objects to be rendered
std::vector<double> load_plane(double a, double b, double c, double d) {
    float plane_a = (float)a;
    float plane_b = (float)b;
    float plane_c = (float)c;
    float plane_d = (float)d;
    
    // Create the equation based on a plane
    // Ax + By + Cz + D = 0
    // All combos of x, y, z that satisfy this equation are on the plane


    std::vector<double> scene_obj;
    // Store a number to indicate shape type
    // 0 is sphere, 1 is sun, 2 is plane, others... TODO
    scene_obj.push_back(2);
    // Then the coordinates of the sphere
    
    scene_obj.push_back(plane_a); //x
    scene_obj.push_back(plane_b); //y
    scene_obj.push_back(plane_c); //z
    scene_obj.push_back(plane_d); //radius
    // Then the RGB (TODO and other things like shininess)
    scene_obj.push_back(r); //Current R
    scene_obj.push_back(g); //Current G
    scene_obj.push_back(b); //Current B
    scene_obj.push_back(1); //Current alpha //TODO How to implement alpha for optional things?

    return scene_obj;
}

// Load a triangle vertex
void get_vertex(float x, float y, float z) {
    //float x = stof(line_tokens[1]);
    //float y = stof(line_tokens[2]);
    //float z = stof(line_tokens[3]);

    //std::vector<float> vertex{x, y, z};
    Eigen::Vector3d vertex(x, y, z);

    vertices.push_back(vertex);
    
}

// Add a triangle to the list of objects to be rendered
std::vector<double> load_triangle(std::vector<std::string> line_tokens) {
    //int vert0_i = stoi(line_tokens[1]) >= 0 ? stoi(line_tokens[1]) - 1 : vertices.size() + stoi(line_tokens[1]);
    //int vert1_i = stoi(line_tokens[2]) >= 0 ? stoi(line_tokens[2]) - 1 : vertices.size() + stoi(line_tokens[2]);
    //int vert2_i = stoi(line_tokens[3]) >= 0 ? stoi(line_tokens[3]) - 1 : vertices.size() + stoi(line_tokens[3]);
    
    // Load the vertices
    float x0 = stof(line_tokens[2]);
    float y0 = stof(line_tokens[3]);
    float z0 = stof(line_tokens[4]);
    // Save to vertices list
    get_vertex(x0, y0, z0);
    int vert0_i = vertices.size() - 1;
    
   
    float x1 = stof(line_tokens[5]);
    float y1 = stof(line_tokens[6]);
    float z1 = stof(line_tokens[7]);
    // Save to vertices list
    get_vertex(x1, y1, z1);
    int vert1_i = vertices.size() - 1;

    float x2 = stof(line_tokens[8]);
    float y2 = stof(line_tokens[9]);
    float z2 = stof(line_tokens[10]);  
    // Save to vertices list
    get_vertex(x2, y2, z2);
    int vert2_i = vertices.size() - 1;
    
    // Vertices list may have duplicates

    
    // Create the triangle based on the plane creation method

    Eigen::Vector3d vert0(x0, y0, z0);
    Eigen::Vector3d vert1(x1, y1, z1);
    Eigen::Vector3d vert2(x2, y2, z2);

    Eigen::Vector3d ab = vert0 - vert1;
    Eigen::Vector3d ac = vert0 - vert2;

    Eigen::Vector3d normal = ab.cross(ac);

    

    float plane_a = normal[0];
    float plane_b = normal[1];
    float plane_c = normal[2];
    float plane_d = -1 * (normal[0] * vert2[0]) - normal[1] * vert2[1] - normal[2] * vert2[2];


    std::vector<double> scene_obj;
    // Store a number to indicate shape type
    // 0 is sphere, 1 is sun, 2 is plane, others... TODO
    // 3 is triangle
    scene_obj.push_back(3);
    // Then the coordinates of the sphere
    
    scene_obj.push_back(plane_a); //x
    scene_obj.push_back(plane_b); //y
    scene_obj.push_back(plane_c); //z
    scene_obj.push_back(plane_d); //radius
    // Then the RGB (TODO and other things like shininess)
    scene_obj.push_back(r); //Current R
    scene_obj.push_back(g); //Current G
    scene_obj.push_back(b); //Current B
    scene_obj.push_back(1); //Current alpha //TODO How to implement alpha for optional things?

    // Store the vertices of the triangle
    // This only stores the indices of the vertices, not the actual vertex values
    scene_obj.push_back((double)vert0_i);
    scene_obj.push_back((double)vert1_i);
    scene_obj.push_back((double)vert2_i);

    std::cout << scene_obj[0] << std::endl;
    std::cout << scene_obj[1] << std::endl;
    std::cout << scene_obj[2] << std::endl;
    std::cout << scene_obj[3] << std::endl;
    std::cout << scene_obj[4] << std::endl;
    std::cout << scene_obj[5] << std::endl;
    std::cout << scene_obj[6] << std::endl;

    return scene_obj;
}

// Add a triangle to the list of objects to be rendered
std::vector<double> make_triangle(float x0, float y0, float z0, float x1, float y1, float z1, float x2, float y2, float z2) {
    //int vert0_i = stoi(line_tokens[1]) >= 0 ? stoi(line_tokens[1]) - 1 : vertices.size() + stoi(line_tokens[1]);
    //int vert1_i = stoi(line_tokens[2]) >= 0 ? stoi(line_tokens[2]) - 1 : vertices.size() + stoi(line_tokens[2]);
    //int vert2_i = stoi(line_tokens[3]) >= 0 ? stoi(line_tokens[3]) - 1 : vertices.size() + stoi(line_tokens[3]);

    // Save to vertices list
    get_vertex(x0, y0, z0);
    int vert0_i = vertices.size() - 1;
   
    // Save to vertices list
    get_vertex(x1, y1, z1);
    int vert1_i = vertices.size() - 1;

    // Save to vertices list
    get_vertex(x2, y2, z2);
    int vert2_i = vertices.size() - 1;
    
    // Vertices list may have duplicates
    
    // Create the triangle based on the plane creation method

    Eigen::Vector3d vert0(x0, y0, z0);
    Eigen::Vector3d vert1(x1, y1, z1);
    Eigen::Vector3d vert2(x2, y2, z2);

    Eigen::Vector3d ab = vert0 - vert1;
    Eigen::Vector3d ac = vert0 - vert2;

    Eigen::Vector3d normal = ab.cross(ac);

    

    float plane_a = normal[0];
    float plane_b = normal[1];
    float plane_c = normal[2];
    float plane_d = -1 * (normal[0] * vert2[0]) - normal[1] * vert2[1] - normal[2] * vert2[2];


    std::vector<double> scene_obj;
    // Store a number to indicate shape type
    // 0 is sphere, 1 is sun, 2 is plane, others... TODO
    // 3 is triangle
    scene_obj.push_back(3);
    // Then the coordinates of the sphere
    
    scene_obj.push_back(plane_a); //x
    scene_obj.push_back(plane_b); //y
    scene_obj.push_back(plane_c); //z
    scene_obj.push_back(plane_d); //radius
    // Then the RGB (TODO and other things like shininess)
    scene_obj.push_back(r); //Current R
    scene_obj.push_back(g); //Current G
    scene_obj.push_back(b); //Current B
    scene_obj.push_back(1); //Current alpha //TODO How to implement alpha for optional things?

    // Store the vertices of the triangle
    // This only stores the indices of the vertices, not the actual vertex values
    scene_obj.push_back((double)vert0_i);
    scene_obj.push_back((double)vert1_i);
    scene_obj.push_back((double)vert2_i);


    return scene_obj;
}

// Add a quadrilateral object (not an infinite plane)
std::vector<std::vector<double>> load_quad(std::vector<std::string> line_tokens) {
    // Cut this into 2 triangles
    float x0 = stof(line_tokens[2]);
    float y0 = stof(line_tokens[3]);
    float z0 = stof(line_tokens[4]);
   
    float x1 = stof(line_tokens[5]);
    float y1 = stof(line_tokens[6]);
    float z1 = stof(line_tokens[7]);

    float x2 = stof(line_tokens[8]);
    float y2 = stof(line_tokens[9]);
    float z2 = stof(line_tokens[10]);

    float x3 = stof(line_tokens[11]);
    float y3 = stof(line_tokens[12]);
    float z3 = stof(line_tokens[13]);
    
    // Find the diagonals and pick one as the long edge to "triangulate" on
    // Get distance between the points
    Eigen::Vector3d vertA(x0, y0, z0);
    Eigen::Vector3d vertB(x1, y1, z1);
    Eigen::Vector3d vertC(x2, y2, z2);
    Eigen::Vector3d vertD(x3, y3, z3);
    
    

    // Declare a holder for all triangle objects we will generate
    std::vector<std::vector<double>> triangle_objs;
    
    // First triangle is (x0,y0,z0), (x1,y1,z1), (x3,y3,z3)
    triangle_objs.push_back(make_triangle(x0, y0, z0, x1, y1, z1, x2, y2, z2));
    // Second triangle is (x1,y1,z1), (x2,y2,z2), (x3,y3,z3) 
    //triangle_objs.push_back(make_triangle(x1, y1, z1, x2, y2, z2, x3, y3, z3));
    triangle_objs.push_back(make_triangle(x2, y2, z2, x3, y3, z3, x0, y0, z0));
    triangle_objs.push_back(make_triangle(x3, y3, z3, x0, y0, z0, x1, y1, z1));
    
    // Feed this into load_plane
    return triangle_objs;
}



void process_scene(std::vector<unsigned char> &image, unsigned width, unsigned height, std::vector<std::vector<double>> scene_objs) {
    // Loop over every pixel to see if it intersects an object
    // Take the first object it intersects as the one driving the pixel's new output

    double sx;
    double sy;

    // Declare the right
    // Default is (1,0,0)
    Eigen::Vector3d right(1, 0, 0);

    // Goal: Make the movable vector, up, perpendicular to the known vector, forward
    // 1: Find right so that it's perpendicular to both forward and up
    right = forward.cross(up);
    right.normalize();
    
    // 2: Change up to be perpendicular to forward and right
    up = right.cross(forward); 
    up.normalize();

    // Get image area
    unsigned max_w_h = std::max(width, height);

    float x_f;
    float y_f;

    for (unsigned y = 0; y < height; y++) {
        for (unsigned x = 0; x < width; x++) {
            // Begin aa loop
            Eigen::Vector3d final_rgb_avg(0,0,0);
            
            std::vector<Eigen::Vector4d> all_rgbs;
            
            double r_fish;
            Eigen::Vector3d forward_fish;
            
            // Background color
            // Will only show if hit by a ray
            Eigen::Vector4d final_rgb_alpha(0, 0, 0, 0);

            x_f = x;
            y_f = y;
            
            // Calculate sx
            sx = (2 * (double)x_f - width) / (double)(max_w_h);

            // Calculate sy
            sy = (height - (2 * (double)y_f)) / (double)(max_w_h);

            // Find distance from at to eye
            //Eigen::Vector3d dist_vec = at - eye;
            double dist = 1;// dist_at_to_eye;
            
            // Calculate viewlength
            double viewheight = 2 * (dist) * std::tan(fov_rad/2);
            
            double viewwidth = viewheight * (width/height);
            
            // Divide sx and sy by viewlength
            //sx = sx / viewlength;
            //sy = sy / viewlength;
            sx = sx * viewwidth;
            sy = sy * viewheight;
            forward_fish = forward;
            r_fish = sqrt(pow(sx, 2) + pow(sy, 2));
             
            // Do not shoot rays if r_fish > 1
            if (r_fish <= 1) {
                // Break the origin and direction of this pixel down
                // Origin
                Eigen::Vector3d origin_v = eye;
                Eigen::Vector3d direction_v;
                direction_v = (sqrt(1 - pow(r_fish, 2)) * forward_fish) + (sx * right) + (sy * up);

                std::vector<std::vector<double>> solutions = get_intersect_all(origin_v, direction_v, scene_objs);
                std::vector<double> solution = get_min_solution(solutions);

                int i_min = solution[0];
                double t_min = solution[1];

                // Now that we've finished looping through all objects that this pixel's
                // // ray may intersect, we can see what to color this pixel
                // Do not color an object if the intersection is "behind" the origin
                if (i_min >= 0) { // An object was found to intersect with this pixel's ray 
                // TODO Better way to avoid drawing intersections "behind" the ray's origin
                // // What does "behind" mean in this context?

                    // Color the pixel based on this object's RGB value
                    // Get the diffuse color of this object
                    float this_r = scene_objs[i_min][5]; // R
                    float this_g = scene_objs[i_min][6]; // G
                    float this_b = scene_objs[i_min][7]; // B
                    float this_a = scene_objs[i_min][8]; // A
                    
                    std::cout << "This object color is " << this_r << ", " << this_g << ", " << this_b << std::endl;

                    // Turn the rgb into vector (don't include alpha)
                    Eigen::Vector3d obj_color = Eigen::Vector3d(this_r, this_g, this_b);

                    // The color will be the sum of, over all lights:
                    // // (object color) * (light color) * (normal (dot) (direction to light))
                    // First, find the normal at this intersection
                    // // This normal vector is perpendicular to the surface
                    // // at the point where this pixel's ray hit it.
                    // Use t_min to find the x,y,z in space where it hit.
                    Eigen::Vector3d collision_pt = origin_v + (t_min * direction_v);

                    Eigen::Vector3d n;
                    // Calculate the normals
                    // Normals will be different depending on the kind of object this is.
                    if (scene_objs[i_min][0] == 0) { // This is a sphere
                        // The normal will be the collision_pt - center (normal is pointing straight out of sphere center)
                        n = collision_pt - Eigen::Vector3d(scene_objs[i_min][1], scene_objs[i_min][2], scene_objs[i_min][3]);
                    } else if (scene_objs[i_min][0] == 2 || scene_objs[i_min][0] == 3) { // This is a plane or triangle
                        // Use barycentric coordinates to interpolate the normals.
                        if (scene_objs[i_min][0] == 3 && read_normals) {
                            // Get the vertices of the triangle
                            Eigen::Vector3d vert0 = vertices[(int)scene_objs[i_min][9]];
                            Eigen::Vector3d vert1 = vertices[(int)scene_objs[i_min][10]];
                            Eigen::Vector3d vert2 = vertices[(int)scene_objs[i_min][11]];

                            // Get the normals of each vertex
                            Eigen::Vector3d norm0 = normals[(int)scene_objs[i_min][9]];
                            Eigen::Vector3d norm1 = normals[(int)scene_objs[i_min][10]];
                            Eigen::Vector3d norm2 = normals[(int)scene_objs[i_min][11]];

                            // Turn this point into a "vector"
                            Eigen::Vector3d this_point = collision_pt;

                            // Get vectors perpendicular to the edges
                            Eigen::Vector3d ab = vert0 - vert1;
                            Eigen::Vector3d ac = vert0 - vert2;
                            Eigen::Vector3d bc = vert1 - vert2;

                            // Find the normal
                            // Either of the below works
                            //Eigen::Vector3d norm_bary = ab.cross(ac);
                            Eigen::Vector3d norm_bary(scene_objs[i_min][1], scene_objs[i_min][2], scene_objs[i_min][3]);

                            Eigen::Vector3d a1 = (vert2 - vert0).cross(norm_bary);
                            Eigen::Vector3d a2 = (vert1 - vert0).cross(norm_bary);

                            Eigen::Vector3d e1 = a1 / a1.dot(vert1 - vert0);
                            Eigen::Vector3d e2 = a2 / a2.dot(vert2 - vert0);

                            double w2 = e2.dot(this_point - vert0);
                            double w1 = e1.dot(this_point - vert0);
                            double w0 = 1 - w1 - w2;
                            
                            // After getting final barycentric coordinates, multiply by normal
                            // // to get interpolated normal
                            n = (norm0 * w0) + (norm1 * w1) + (norm2 * w2);
                        } else {
                            n = Eigen::Vector3d(scene_objs[i_min][1], scene_objs[i_min][2], scene_objs[i_min][3]);
                            //std::cout << "Normals: " << n.normalized()(0) << " " << n.normalized()(1) << " " << n.normalized()(2) << std::endl;
                        }
                    }

                    n.normalize();
                    
                    
                    // Invert normal if it points away from eye
                    // This will make objects two-sided
                    std::cout << "Direction_v is " << direction_v << std::endl;
                    if (direction_v.dot(n) > 0) {
                        std::cout << "Normals wrong: " << n(0) << " " << n(1) << " " << n(2) << std::endl;
                        std::cout << "Switch the normal to fix." << std::endl;
                        n = -1 * n;
                    }
                    std::cout << "Normals: " << n(0) << " " << n(1) << " " << n(2) << std::endl;
                    
                    // Make holder for final pixel color
                    Eigen::Vector3d final_rgb(0,0,0);

                    // Grab anything in scene_objs with code == 1 (sun) or 4 (bulb)
                    for (int i = 0; i < scene_objs.size(); i++) {
                        if (scene_objs[i][0] == 1) { // This is a sun object
                            Eigen::Vector3d sun_dir = Eigen::Vector3d(scene_objs[i][1], scene_objs[i][2], scene_objs[i][3]);
                            // We need to normalize this light direction
                            sun_dir.normalize();
                            Eigen::Vector3d sun_color = Eigen::Vector3d(scene_objs[i][4], scene_objs[i][5], scene_objs[i][6]);

                            // Shoot a shadow ray from the intersection point to the light source
                            // If the ray intersects something before reaching the light source, this will be in shadow
                            // Else, color is influenced by light's color
                            Eigen::Vector3d new_collision_pt = move_origin_for_shadow(collision_pt, sun_dir);

                            std::vector<std::vector<double>> solutions_shadow = get_intersect_all(new_collision_pt, sun_dir, scene_objs);

                            std::vector<double> solution_shadow = get_min_solution_shadow(solutions_shadow, i_min);
                            
                            if ((int)solution_shadow[0] >= 0) { // There was an intersection
                                final_rgb = final_rgb; //TODO Don't do anything if this sun is casting a shadow
                                //final_rgb = (n * 0.5);
                                //final_rgb.array() += 0.5;
                                //double sun_dot = n.dot(sun_dir);
                                //final_rgb = final_rgb + obj_color.cwiseProduct(sun_color * std::max((double)0, sun_dot));
                            } else { //No intersection found
                                double sun_dot = n.dot(sun_dir);
                                final_rgb = final_rgb + obj_color.cwiseProduct(sun_color * std::max((double)0, sun_dot));
                                //final_rgb = (n * 0.5);
                                //final_rgb.array() += 0.5;
                            }

                        } else if (scene_objs[i][0] == 4) { // This is a bulb object
                            // Get position of bulb in 3D space
                            Eigen::Vector3d bulb_pos = Eigen::Vector3d(scene_objs[i][1], scene_objs[i][2], scene_objs[i][3]);
                            Eigen::Vector3d bulb_dir = bulb_pos - collision_pt;
                            double dist = bulb_dir.norm();
                            bulb_dir.normalize();
                            Eigen::Vector3d bulb_color = Eigen::Vector3d(scene_objs[i][4], scene_objs[i][5], scene_objs[i][6]);
                            
                            // Shoot a shadow ray from the intersection point to the light source
                            // If the ray intersects something before reaching the light source, this will be in shadow
                            // Else, color is influenced by light's color
                            Eigen::Vector3d new_collision_pt = move_origin_for_shadow(collision_pt, bulb_dir);      

                            std::vector<std::vector<double>> solutions_shadow = get_intersect_all(new_collision_pt, bulb_dir, scene_objs);

                            std::vector<double> solution_shadow = get_min_solution_shadow(solutions_shadow, i_min);

                            if ((int)solution_shadow[0] >= 0) { // There was an intersection
                                final_rgb = obj_color; //final_rgb; //TODO Don't do anything if there's an object between vector from intersection to bulb, casting a shadow
                                
                            } else { //No intersection found
                                // Find the diffuse amount
                                double bulb_dot = n.dot(bulb_dir);
                                double diffuseAmt = std::max((double)0, bulb_dot);
                                
                                // See effect with distance from light
                                double attenuation = 1 / (dist * dist);
                                
                                //Ambient. Assume Ambient coeffient, ka, is close to 0, so 0.1
                                Eigen::Vector3d ambient = bulb_color * 0.1;
                                
                                // Diffuse (use object's kd)
                                Eigen::Vector3d diffuse = bulb_color * kd * std::max((double)0, bulb_dot);
                                
                                
                                // Specular (use object's ks)
                                // To Do
                                Eigen::Vector3d specular = bulb_color * ks; // To Do
                                
                                // Final pixel color is ambient + diffuse + specular
                                //final_rgb = final_rgb + obj_color;
                                //final_rgb = final_rgb + obj_color.cwiseProduct(bulb_color * std::max((double)0, bulb_dot)) ; //obj_color.cwiseProduct(sun_color * std::max((double)0, bulb_dot));
                                final_rgb = final_rgb + obj_color.cwiseProduct((ambient + diffuse) * attenuation) ;
                                //final_rgb = (n * 0.5);
                                //final_rgb.array() += 0.5;
                            }                            
                            
                            // We need to normalize the direction from this vertex to the light
                            //(bulb_dir - vertex_loc)
                        }
                            
                            
                    }
                    // Save this final_rgb for averaging later
                    Eigen::Vector4d final_rgb_alpha_plus(final_rgb[0], final_rgb[1], final_rgb[2], this_a);
                    
                    final_rgb_alpha = final_rgb_alpha + final_rgb_alpha_plus;
                    std::cout << "final_rgb for this pixel is " << final_rgb_alpha[0] << ", " << final_rgb_alpha[1] << ", " << final_rgb_alpha[2] << std::endl;
                } else { // Ray didn't hit anything, so show background color
                    final_rgb_alpha = bg_c;
                }
                //all_rgbs.push_back(final_rgb_alpha);
                
                unsigned this_r_u = convert_linear_color_to_unsigned_srgb(final_rgb_alpha[0]);
                unsigned this_g_u = convert_linear_color_to_unsigned_srgb(final_rgb_alpha[1]);
                unsigned this_b_u = convert_linear_color_to_unsigned_srgb(final_rgb_alpha[2]);
                unsigned this_a_u = convert_linear_color_to_unsigned_linear(final_rgb_alpha[3]);            

                color_pixel(image, width, x, y, this_r_u, this_g_u, this_b_u, this_a_u);            
            
            }
               
        }
        
    }
}

// From LodePNG CPP example (LodePNG's manual)
// https://raw.githubusercontent.com/lvandeve/lodepng/master/examples/example_encode.cpp
unsigned save_image_png(std::string filename, std::vector<unsigned char> image, unsigned width, unsigned height) {
    std::vector<unsigned char> png;

    unsigned error = lodepng::encode(png, image, width, height);
    if(!error) {
        lodepng::save_file(png, filename);
        return 0; //Return code for success
    } else { // There's an error. Print it.
        if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
        return 1; // Return code for failure
    }

}


int main(int argc, char *argv[]){
    // Open the file into an array of each line
    
    const char* input_filename = argv[1];
    
    // Read in NFF file format
    std::vector<std::string> lines = text_file_to_lines(input_filename);

    // Declare a holder for the scene's objects
    std::vector<std::vector<double>> scene_objs;
    
    for (std::string line : lines) {
        std::vector<std::string> line_tokens = line_to_tokens(line);
        std::cout << line << "\n";
            
        if (line_tokens[0] == "v") {
            // Declaring view parts. Do nothing
        } else if ((line_tokens[0]) == "from") { // DONE Translated from "from"
            get_eye_info(line_tokens);
        } else if ((line_tokens[0]) == "at") { // DONE
            get_at_info(line_tokens);
        } else if ((line_tokens[0]) == "up") { // DONE
            get_up_info(line_tokens);
        } else if ((line_tokens[0]) == "angle") { // DONE
            get_angle_info(line_tokens);
        } else if ((line_tokens[0]) == "hither") { // To Do, incorporate in scene draw code
            get_hither_info(line_tokens);
        } else if ((line_tokens[0]) == "resolution") { // DONE
            get_resolution_info(line_tokens);
        } else if ((line_tokens[0]) == "l") { // DONE Translated from "sun" but fixed to "bulb"
            scene_objs.push_back(load_bulb(line_tokens));
        } else if ((line_tokens[0]) == "s") { // DONE Translated from "sphere"
            scene_objs.push_back(load_sphere(line_tokens));
        } else if ((line_tokens[0]) == "b") { // DONE
            get_bg_info(line_tokens);
        } else if ((line_tokens[0]) == "f") { // DONE, translated from "color" //TODO Implement non-RGB parts
            color(r, g, b, line_tokens);
        } else if ((line_tokens[0]) == "p" && line_tokens[1] == "3") { // To do
            scene_objs.push_back(load_triangle(line_tokens));
        } else if ((line_tokens[0]) == "p" && line_tokens[1] == "4") { // To do
            // Export a list of triangles
            // For each member of the list, push back the triangle object to scene_objs
            std::vector<std::vector<double>> quad_tris = load_quad(line_tokens);
            for (std::vector<double> tri : quad_tris) {
                scene_objs.push_back(tri);
            }
            
        }
        
    }

    // Instantiate the image
    std::vector<unsigned char> image;
    image.resize(width * height * 4);

    // After reading in all commands, process how the rays interact with the 
    // // objects loaded into scene_objs
    // For each pixel, see if its ray intersects with any objects
    // If so, make the pixel take on the color of the object.
    process_scene(image, width, height, scene_objs);

    // Declare output file name
    std::string output_filename = "raytracer_output.png";
    
    // Save the image
    save_image_png(output_filename, image, width, height);
    std::cout << "Saved to " << output_filename << "\n";

    return 0;
    
} 