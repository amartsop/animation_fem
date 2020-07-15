
void demo_animation() {

    Gnuplot gp;

    std::cout << "Press Ctrl-C to quit (closing gnuplot window doesn't quit)." << std::endl;


    
    float fs = 1;
    
    auto start = std::chrono::steady_clock::now();
    
    float t = 0;

    gp << "set style line 1 lc rgb 'black' lw 3 pt 7 ps 1.5\n";

    gp << "set yrange [-5:5]\n";
    gp << "set xrange [-5:5]\n";
    gp << "set zrange [-5:5]\n";
    // gp << "set hidden3d\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'z'\n";
    gp << "set view 60, 45, 1, 1\n";
    gp << "set grid\n";

    while(1)
    {
   
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<float> elapsed_seconds = end-start;

    // Real time
    float t = t + elapsed_seconds.count();

    // Model
    float x1 = 0;
    float x2 = 0;
    float x3 = 0;
    
    float y1 = -1;
    float y2 = 0;
    float y3 = 1;
    
    nfloat z1 = 0;
    float z2 = sinf(2*M_PIf32 * fs * t);
    float z3 = sinf(2*M_PIf32 * fs * t);
    
    arma::fmat xy = { {x1, y1, z1}, {x2, y2, z2}, {x3, y3, z3} };


    gp << "splot '-' with linespoints ls 1\n";
    gp.send1d(xy);
    gp.flush();

    // Update time
    start = end;

    }
    
    
}

void demo_cylinder(void)
{

    
    arma::fvec roa_fe_fe = {0, 0, 0};
    arma::fvec theta_vec = {0, 0, 0};

    float radius = 1;
    float height = 3;
    int side_count = 10;

    arma::fmat vertices = arma::zeros<arma::fmat>(2 * side_count, 3);

    // Cylinder vertices
    for (int i = 0; i < side_count; i++)
    {
        float theta = 2 * (M_PIf32 / side_count) * i;
        arma::fvec vert_row1 = {radius * cosf(theta), radius * sinf(theta), 0}; 
        arma::fvec vert_row2 = {radius * cosf(theta), radius * sinf(theta), height}; 
        vertices.row(i) = vert_row1.t();
        vertices.row(side_count + i) = vert_row2.t();
    }

    // Cylinder side faces
    arma::imat side_faces = arma::zeros<arma::imat>(side_count, 4);
    for (int i = 0; i < side_count - 1; i++)
    {
        arma::ivec side_faces_row = {i + 1, i + 2, side_count + i + 2, 
            side_count + i + 1}; 
        side_faces.row(i) = side_faces_row.t();
    }
    arma::ivec side_faces_last_row = {side_count, 1, side_count + 1, 2*side_count};
    side_faces.row(side_count - 1) = side_faces_last_row.t();

    // Cylinder Bottom faces
    arma::ivec lower_bottom_face = arma::regspace<arma::ivec>(1, side_count);
    arma::ivec upper_bottom_face = arma::regspace<arma::ivec>(side_count + 1,  
        2 * side_count);

    arma::imat bottom_faces = arma::zeros<arma::imat>(2, side_count);
    bottom_faces.row(0) = lower_bottom_face.t();
    bottom_faces.row(1) = upper_bottom_face.t();


    // DRAWING
    arma::fmat ad = {{0,0,0}, {1,0,0}, {0.5,2,0}};

    std::string x0 = std::to_string(ad(0, 0));
    std::string y0 = std::to_string(ad(0, 1));
    std::string z0 = std::to_string(ad(0, 2));
    std::string coord0 = x0 + "," + y0 + "," + z0;

    std::string x1 = std::to_string(ad(1, 0));
    std::string y1 = std::to_string(ad(1, 1));
    std::string z1 = std::to_string(ad(1, 2));
    std::string coord1 = x1 + "," + y1 + "," + z1;
    
    std::string x2 = std::to_string(ad(2, 0));
    std::string y2 = std::to_string(ad(2, 1));
    std::string z2 = std::to_string(ad(2, 2));
    std::string coord2 = x2 + "," + y2 + "," + z2;

    Gnuplot gp;

    // gp << "set style line 1 lc rgb 'black' lw 3 pt 7 ps 1.5\n";

    gp << "set yrange [-5:5]\n";
    gp << "set xrange [-5:5]\n";
    gp << "set zrange [-5:5]\n";
    gp << "set hidden3d\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'z'\n";
    gp << "set view 60, 45, 1, 1\n";
    gp << "set grid\n";




    gp << "set object 1 polygon from " << coord0 << " to " << coord1 << " to " <<
        coord2 << " to " << coord0 << " fillstyle transparent solid 0.5\n";

    gp << "splot 1\n";

    // gp.flush();


    
    // Draw Cylinder
    // gp << "set dgrid3d\n";
    // gp << "set hidden3d\n";
    // gp << "splot '-' with lines\n";
    // gp.send1d(ad);
    // gp.flush();

    // gp << "splot '-' with filledcurves\n";
    // gp.send1d(ad);
    // gp.flush();


}