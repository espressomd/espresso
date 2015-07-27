class Geometry {
  public:
    void set_box_l(Vector3d x) {m_box_l =x;}
    void set_my_left(Vector3d x) {m_my_left=x;}
    void set_my_right(Vector3d x) { m_my_right=x;}
    Vector3d get_box_l() {return m_box_l;}
    Vector3d get_my_left() { return m_my_left; }
    Vector3d get_my_right() { return m_my_right; }
  private:
    Vector3d m_box_l;
    Vector3d m_my_left;
    Vector3d m_my_right;
};

extern Geometry* geometry;


