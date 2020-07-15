#ifndef NEEDLE_GEOMETRY_H
#define NEEDLE_GEOMETRY_H

class NeedleGeometry
{
public:
    static float get_handle_radius(void) { return m_radius; }
    static float get_handle_height(void) { return m_height; }
    static int get_handle_side_count(void) { return m_side_count; }

    static float get_needle_offset_z(void) {return m_lz; }
    static float get_needle_offset_x(void) {return m_lx; }
    static float get_needle_length(void) { return m_ln; }

protected:
    // Handle geometry
    static constexpr float m_radius = 0.017; // Handle radius (m) 
    static constexpr float m_height = 0.153; // Handle height (m)
    static constexpr int m_side_count = 20; // Side numbers for cylinder
    
    // Needle geometry
    static constexpr float m_lz = 0.00767; // Needle offset z direction (m)
    static constexpr float m_lx = 0.0765; // Needle offset x direction (m)
    static constexpr float m_ln = 0.207; // Needle length (m)
};




#endif