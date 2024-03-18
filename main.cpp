#define _USE_MATH_DEFINES
#include <chrono>
#include "M_matrices.h"
#include "Color_picker.h"
#include <Windows.h>
#include <conio.h>
#include <string>

using namespace std;

struct Point;
struct Arrete;
struct Object;
struct Environnement;
struct Camera;
Arrete point_to_arrrete(Point p);
Point arrete_to_point(Arrete a);
bool key_pressed(Environnement &E);

const size_t LIGNES(41), COL(101);

struct Point
{
    Point(){}
    Point(double x, double y, double z): m_x(x), m_y(y), m_z(z){}
    //Point(point): m_x(point.m_x), m_y(point.m_y), m_z(point.m_z){}
    double m_x;
    double m_y;
    double m_z;

    double norme()
    {
        return sqrt(pow(this->m_x,2)+pow(this->m_y,2)+pow(this->m_z,2));
    }
    void set_coord(double x, double y, double z)
    {
        m_x = x;
        m_y = y;
        m_z = z;
    }
    void set_coord(vector<double> coord)
    {
        m_x = coord[0];
        m_y = coord[1];
        m_z = coord[2];
    }
    M_matrices to_matrice()
    {
        vector<double> mat_vect({m_x, m_y, m_z});
        M_matrices mat(3,1,mat_vect);
        return mat;
    }
    void operator=(vector<double> coord)
    {
        m_x = coord[0];
        m_y = coord[1];
        m_z = coord[2];
    }
    void operator+=(Point p)
    {
        m_x += p.m_x;
        m_y += p.m_y;
        m_z += p.m_z;
    }
    vector<double> point_to_vector()
    {
        vector<double> v;
        v.push_back(m_x);
        v.push_back(m_y);
        v.push_back(m_z);
        return v;
    }
    double distance_to_point(Point B)
    {
        return sqrt(pow(B.m_x-this->m_x,2)+pow(B.m_y-this->m_y,2)+pow(B.m_z-this->m_z,2));
    }
};

Point matrice_to_point(M_matrices mat)
{
    if (mat.get_matrice_vect().size() == 3)
        return Point(mat.get_matrice_vect()[0], mat.get_matrice_vect()[1], mat.get_matrice_vect()[2]);
    else
        cerr<<"matrice_to_point : size_error"<<endl;
        return Point(0,0,0);
}

Point operator*(Point a, double x)
{
    Point S;
    S.m_x = a.m_x *x;
    S.m_y = a.m_y *x;
    S.m_z = a.m_z *x;
    return S;
}

Point operator/(Point a, double x)
{
    Point S;
    S.m_x = a.m_x /x;
    S.m_y = a.m_y /x;
    S.m_z = a.m_z /x;
    return S;
}

Point operator+(Point a, Point b)
{
    Point S;
    S.m_x = a.m_x + b.m_x;
    S.m_y = a.m_y + b.m_y;
    S.m_z = a.m_z + b.m_z;
    return S;
}

Point operator-(Point a, Point b)
{
    Point S;
    S.m_x = a.m_x - b.m_x;
    S.m_y = a.m_y - b.m_y;
    S.m_z = a.m_z - b.m_z;
    return S;
}



struct Arrete
{
    Arrete(){}
    Arrete(double x, double y, double z): m_x(x), m_y(y), m_z(z){}
    double m_x;
    double m_y;
    double m_z;

    double get_norm()
    {
        return sqrt(pow(m_x,2)+pow(m_y,2)+pow(m_z,2));
    }

    void transform_to_unit_vect()
    {
        double norme = get_norm();
        m_x /= norme;
        m_y /= norme;
        m_z /= norme;
    }
};

Arrete operator+(Arrete a, Arrete b)
{
    Arrete c;
    c.m_x = a.m_x + b.m_x;
    c.m_y = a.m_y + b.m_y;
    c.m_z = a.m_z + b.m_z;
    return c;
}
Arrete operator*(Arrete a, double b)
{
    Arrete c;
    c.m_x = a.m_x * b;
    c.m_y = a.m_y * b;
    c.m_z = a.m_z * b;
    return c;
}
Arrete operator/(Arrete a, double b)
{
    Arrete c;
    c.m_x = a.m_x / b;
    c.m_y = a.m_y / b;
    c.m_z = a.m_z / b;
    return c;
}


struct Object
{
    Object(): m_rot_x(0), m_rot_y(0), m_rot_z(0), m_origine(Point(0,0,0)), m_definig_points({}), m_arretes({}){}
    Object(Point normale_plan, Point position) : m_name("Camera"), m_rot_x(0), m_rot_y(0), m_rot_z(0),
                                    m_origine(Point(0,0,0)), m_direction({normale_plan}), m_arretes({}),
                                    m_definig_points({position, Point(1,0,0), Point(0,1,0), Point(0,0,-1)}){}
    double m_rot_x;
    double m_rot_y;
    double m_rot_z;
    double m_length;
    double m_rayon;
    double m_pas;
    M_matrices m_matrice_rot = get_matrice_identite(3,3);
    Point m_origine;
    vector<Point> m_definig_points;
    Point m_direction;
    vector<Arrete> m_arretes;
    string m_name;

    void make_line(double length, double pas, string name)
    {
        m_name = name;
        m_definig_points.clear();
        m_length = length;
        m_pas = pas;
        Point p = Point(0,0,0);

        for(int i = 0; i*pas<length; i++)
        {
            m_definig_points.push_back(Point(0,0,i*pas));
        }
        m_definig_points.push_back(Point(0,0,length));
        set_centered();
    }
    void make_sphere(double rayon, double pas, string name)
    {
        m_name = name;
        m_definig_points.clear();
        m_rayon = rayon;
        m_length = 2*rayon;//pour l'affichage
        m_definig_points.push_back(Point(rayon,rayon,0));
        m_definig_points.push_back(Point(rayon,-rayon,0));

        for(int i = 0; i*pas<rayon; i++)//de 0 à 2*rayon
        {
            //cout<<"mdr : "<<i<<endl;
            double x_pos = sqrt(pow(rayon,2)-pow(i*pas,2))+rayon;
            double x_neg = -sqrt(pow(rayon,2)-pow(i*pas,2))+rayon;
            m_definig_points.push_back(Point(x_pos,i*pas,0));
            m_definig_points.push_back(Point(x_neg,i*pas,0));

            m_definig_points.push_back(Point(x_pos,-i*pas,0));
            m_definig_points.push_back(Point(x_neg,-i*pas,0));

            m_definig_points.push_back(Point(x_pos,0,i*pas));
            m_definig_points.push_back(Point(x_neg,0,i*pas));

            m_definig_points.push_back(Point(x_pos,0,-i*pas));
            m_definig_points.push_back(Point(x_neg,0,-i*pas));

            m_definig_points.push_back(Point(rayon,i*pas,x_pos-rayon));
            m_definig_points.push_back(Point(rayon,i*pas,x_neg-rayon));

            m_definig_points.push_back(Point(rayon,-i*pas,x_pos-rayon));
            m_definig_points.push_back(Point(rayon,-i*pas,x_neg-rayon));
        }
        set_centered();
    }

    void make_cube(double length, double pas, string name)
    {
        m_name = name;
        m_definig_points.clear();
        m_arretes.clear();
        m_length = length;
        m_pas = pas;
        Arrete a = Arrete(0,0,0);
        Point p = Point(0,0,0);
        rotate_point(p);
        m_definig_points.push_back(p);
        p.m_z += length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        p.m_y += length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        p.m_z-=length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        p.m_x -= length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        p.m_y-= length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        p.m_z += length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        p.m_y += length;
        rotate_point(p);
        a = point_to_arrrete(p - m_definig_points[m_definig_points.size()-1]);
        m_definig_points.push_back(p);
        m_arretes.push_back(a);
        //dernières arrètes
        a = point_to_arrrete(m_definig_points[3]);//1 à 4
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[0];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[5]);//1 à 6
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[0];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[6] - m_definig_points[1]);//2 à 7
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[1];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[7] - m_definig_points[2]);//3 à 8
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[2];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[7] - m_definig_points[4]);//5 à 8
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[4];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }

        for(size_t i(0); i<m_arretes.size()-5; i++)
        {
            for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[i];
                m_arretes[i].transform_to_unit_vect();
                Arrete v = m_arretes[i] *(j*pas);///peut manquer une arrete

                m_definig_points.push_back(o+arrete_to_point(v));
                //cout<<i<<endl;
            }
        }
        //cout<<"ok"<<endl;

        set_centered();
        //cout<<"ok2"<<endl;
    }
    void set_centered()
    {
        set_center();
        Point temp =  Point(0,0,0) - m_origine;
        translate(temp.point_to_vector());
    }
    void set_center()
    {
        for(size_t i(0); i<m_definig_points.size(); i++)
        {
            m_origine.m_x += m_definig_points[i].m_x;
            m_origine.m_y += m_definig_points[i].m_y;
            m_origine.m_z += m_definig_points[i].m_z;
        }
        m_origine.m_x /= m_definig_points.size();
        m_origine.m_y /= m_definig_points.size();
        m_origine.m_z /= m_definig_points.size();
    }
    void refresh_rotation()
    {
        for(int i(0); i<m_definig_points.size(); i++)
        {
            rotate_point(m_definig_points[i]);
            //cout<<i<<" : "<<m_definig_points[i].m_x<<","<<m_definig_points[i].m_y<<","<<m_definig_points[i].m_z<<endl;
        }
        rotate_point(m_direction);
        set_center();
    }
    void translate_point(Point &p, Point translation)
    {
        p += translation;
    }
    void translate(vector<double> translation)
    {
        if (translation.size() > 3)
            cerr<<"translation impossible car trop de coord."<<endl;
        else
        {
            if (m_name == "Camera")
            {
                Point temp = Point(translation[0], translation[1], translation[2]);
                translate_point(m_definig_points[0], temp);
            }
            else
            {
                for(size_t i(0); i<m_definig_points.size(); i++)
                {
                    Point temp = Point(translation[0], translation[1], translation[2]);
                    translate_point(m_definig_points[i], temp);
                }
            }

        }
    }
    void rotate_x(double angle_rad)
    {
        m_rot_x += angle_rad;
        vector<double> matrice_rot_x = {1,0,0,
                                        0,cos(angle_rad),-sin(angle_rad),
                                        0,sin(angle_rad),cos(angle_rad)};
        m_matrice_rot.set_matrice_vect(matrice_rot_x);
        refresh_rotation();
    }
    void rotate_y(double angle_rad)
    {
        m_rot_y += angle_rad;
        vector<double> matrice_rot_y = {cos(angle_rad),0,sin(angle_rad),
                                        0,1,0,
                                        -sin(angle_rad),0,cos(angle_rad)};
        m_matrice_rot.set_matrice_vect(matrice_rot_y);
        refresh_rotation();
    }

    void rotate_z(double angle_rad)
    {
        m_rot_z += angle_rad;
        vector<double> matrice_rot_z = {cos(angle_rad),-sin(angle_rad),0,
                                        sin(angle_rad),cos(angle_rad) ,0,
                                        0             ,0              ,1};
        m_matrice_rot.set_matrice_vect(matrice_rot_z);
        refresh_rotation();
    }
    void rotate_point(Point &p)
    {
        M_matrices m_point(3,1);
        m_point.set_matrice_vect({p.m_x,p.m_y,p.m_z});
        M_matrices new_point = m_matrice_rot*m_point;
        p.set_coord(new_point.get_matrice_vect());
    }

};

struct Camera : Object
{
    //constructeurs
    //Camera():m_vect_normal(Point(0,0,1)), m_origine_plan(Point()), m_size_x(-1), m_size_y(-1){}
    //Camera(double Ox,double Oy,double Oz, Point origine_plan) : m_vect_normal({Ox,Oy,Oz}),
           // m_origine_plan(origine_plan), m_size_x(-1), m_size_y(-1){}
    Camera(Point plan, Point origine_plan) : Object(plan, origine_plan), m_origine_plan(origine_plan),
            m_vect_normal(plan), m_size_x(-1), m_size_y(-1), camera_point({10,0,0}){}
    //attributs
    Point m_vect_normal;//vecteur normal
    Point m_origine_plan;
    double m_size_x;
    double m_size_y;
    Point camera_point;
    //m_definig_points.push_back(vect_e1);
    //this.m_definig_points.push_back(vect_e2);

    //méthodes

    Point transformation_plan(Point point_R3)
    {
        m_vect_normal = m_direction;
        m_origine_plan = m_definig_points[0];
        //cout<<m_definig_points[0].m_x<<","<< m_definig_points[0].m_y<<","<< m_definig_points[0].m_z<<endl;
        camera_point = m_definig_points[0]-((m_vect_normal/m_vect_normal.norme())*10);//changé le - par un +
        double d = -(m_vect_normal.m_x*m_origine_plan.m_x) - m_vect_normal.m_y*m_origine_plan.m_y - m_vect_normal.m_z*m_origine_plan.m_z;
        Point point_plan_ortho_A;
        Point point_plan_ortho_B;
        double k = (m_vect_normal.m_x*point_R3.m_x+m_vect_normal.m_y*point_R3.m_y+m_vect_normal.m_z*point_R3.m_z+d)
            /(pow(m_vect_normal.m_x,2)+pow(m_vect_normal.m_y,2)+pow(m_vect_normal.m_z,2));
        double k_2 = (m_vect_normal.m_x*camera_point.m_x+m_vect_normal.m_y*camera_point.m_y+m_vect_normal.m_z*camera_point.m_z+d)
            /(pow(m_vect_normal.m_x,2)+pow(m_vect_normal.m_y,2)+pow(m_vect_normal.m_z,2));

        //Point A = caméra, point B = objet
        Point vect_AB = point_R3-camera_point;
        double dist_AB = vect_AB.norme();
        point_plan_ortho_B.set_coord(point_R3.m_x-k*m_vect_normal.m_x, point_R3.m_y-k*m_vect_normal.m_y, point_R3.m_z-k*m_vect_normal.m_z);
        point_plan_ortho_A.set_coord(camera_point.m_x-k_2*m_vect_normal.m_x, camera_point.m_y-k_2*m_vect_normal.m_y, camera_point.m_z-k_2*m_vect_normal.m_z);
        double dist_A = point_plan_ortho_A.distance_to_point(camera_point);
        double dist_B = point_plan_ortho_B.distance_to_point(point_R3);
        Point vect_Aprime_Bprime = point_plan_ortho_B-point_plan_ortho_A;
        Point vect_unit_Aprime_Bprime = vect_Aprime_Bprime/vect_Aprime_Bprime.norme();
        double dist_Aprime_inter = sqrt(pow(dist_AB,2)-pow(dist_A+dist_B,2))-sqrt(pow(dist_A*dist_AB/(dist_A+dist_B)-dist_AB,2)-pow(dist_A,2));



        //Point A_B = point_plan_ortho_B - point_plan_ortho_A;
        //double temp = (dist_A+dist_B)/dist_A;
        //double dist_AX = A_B.norme()/temp;
        //Point unit_A_B = A_B/A_B.norme();
        Point point_plan = point_plan_ortho_B+(vect_unit_Aprime_Bprime*dist_Aprime_inter);
        //pos_point_R3 - m*vect_directeur = pos_point_r2 // m>0
        //double m = -((point_plan.m_x - point_R3.m_x) / m_vect_normal.m_x);
        ///TROUVER EQUATION POUR M < 0 => point derrière caméra
        /*Point vect_dist = point_plan-m_origine_plan;
        double delta_dist = sqrt(pow(vect_dist.m_x,2)+pow(vect_dist.m_y,2)+pow(vect_dist.m_z,2));
        M_matrices matrice_pass(3,3);
        vector<double> vect_mat({m_definig_points[1].m_x,m_definig_points[1].m_y,m_definig_points[1].m_z,
                                 m_definig_points[2].m_x,m_definig_points[2].m_y,m_definig_points[2].m_z,
                                 m_vect_normal.m_x, m_vect_normal.m_y, m_vect_normal.m_z});
        matrice_pass.set_matrice_vect(vect_mat);
        M_matrices mat_point(3,1,{point_plan.m_x,point_plan.m_y,point_plan.m_z});
        M_matrices mat_point_plan = matrice_pass.get_inverse()*mat_point;
        vector<double> temp = mat_point_plan.get_matrice_vect();
        point_plan.set_coord(temp);*/

        double x = point_plan.m_x;
        double y = point_plan.m_y;
        double z = point_plan.m_z;
        double u1 = m_definig_points[1].m_x;
        double u2 = m_definig_points[1].m_y;
        double u3 = m_definig_points[1].m_z;
        double v1 = m_definig_points[2].m_x;
        double v2 = m_definig_points[2].m_y;
        double v3 = m_definig_points[2].m_z;
        double w1 = m_definig_points[3].m_x;
        double w2 = m_definig_points[3].m_y;
        double w3 = m_definig_points[3].m_z;

        vector<double> P_mat_vect({u1,v1,w1,u2,v2,w2,u3,v3,w3});
        M_matrices P_mat(3,3,P_mat_vect);

        M_matrices P_mat_inv(P_mat.get_inverse());
        M_matrices mat_point_plan = P_mat_inv*point_plan.to_matrice();
        Point point_cam_plan(matrice_to_point(mat_point_plan));
        //double coeff_a =
        //double coeff_b = (y*u1-u2*x)/(v2*u1-u2*v1);
        //point_plan.set_coord(0, coeff_a, coeff_b);
        /*double x1=0.;
        double x2=1.;
        double x3=0.;
        double y1=0.;
        double y2=0.;
        double y3=1.;
        double d1=m_definig_points[0].distance_to_point(point_plan);
        double d2=m_definig_points[1].distance_to_point(point_plan);
        double d3=m_definig_points[2].distance_to_point(point_plan);
        vector<double> vect_mat_a = {2*(x3-x1),2*(y3-y1),2*(x3-x2),2*(y3-y2)};
        double temp1 = pow(d1,2)-pow(d3,2)+pow(x3,2)-pow(x1,2)+pow(y3,2)-pow(y1,2);
        double temp2 = pow(d2,2)-pow(d3,2)+pow(x3,2)-pow(x2,2)+pow(y3,2)-pow(y2,2);
        vector<double> vect_mat_b = {temp1,temp2};
        M_matrices mat_a(2,2,vect_mat_a);
        M_matrices mat_b(2,1,vect_mat_b);
        M_matrices mat_res(2,1);
        mat_res = mat_a.get_inverse()*mat_b;
        vector<double> temp = mat_res.get_matrice_vect();
        point_plan.set_coord(temp);*/
        /*if(m < 0 || delta_dist>200)//derrière la caméra ou hors du champ
        {
            return Point(0,0,0);
        }
        else
        {
            return point_plan;
        }*/
        if(camera_point.distance_to_point(point_R3)>m_definig_points[0].distance_to_point(point_R3))
        {
            return point_cam_plan;
        }
        else
        {
            return Point(1000,1000,1000);///a rechanger
        }

    }

};

struct Environnement
{
    //Environnement():m_x_axis(0), m_y_axis(0), m_z_axis(0), m_cam(){}
    Environnement(double x_axis, double y_axis, double z_axis): m_x_axis(x_axis), m_y_axis(y_axis), m_z_axis(z_axis),
                    m_cam(Point(-1,0,0),Point(0,0,0)){}
    //static const size_t LIGNES=40, COL=50;
    char display_tab[LIGNES][COL];
    int display_tab_color[LIGNES][COL];
    int m_current_object_index = 0;

    bool is_camera_view = false;
    double m_x_axis;
    double m_y_axis;
    double m_z_axis;
    vector<Object*> list_objects;
    vector<int> list_objects_color;
    vector<Point> list_objects_coord;
    Camera m_cam;
    //vector<string> COLor_pannel = {"bleu","rouge"};

    void select_next_object()
    {
        if (m_current_object_index == list_objects.size()-1)
            m_current_object_index = 0;
        else
            m_current_object_index++;
    }

    void move_object(Point target)//déplace l'objet courant dans l'environnement
    {
        list_objects_coord[m_current_object_index] += target;
    }

    void clear_all()
    {
        for (int i = 0; i<list_objects.size(); i++)
        {
            delete list_objects[i];
        }
        list_objects.clear();
        list_objects_color.clear();
        list_objects_coord.clear();
    }
    void add_object(Object *item, Point coord, string color)
    {
        list_objects.push_back(item);
        int color_int = transform_text_color(color);
        list_objects_color.push_back(color_int);
        list_objects_coord.push_back(coord);
    }
    void set_environnement(double x_axis, double y_axis, double z_axis)
    {
        m_x_axis = x_axis;
        m_y_axis = y_axis;
        m_z_axis = z_axis;
    }
    Camera get_camera()
    {
        return m_cam;
    }
    void clear_display_tab()
    {
        for(size_t i(0); i<LIGNES; i++)
        {
            for(size_t j(0); j<COL; j++)
            {
                display_tab[i][j] = ' ';
                display_tab_color[i][j] = 15;//blanc
            }
        }
    }
    void set_display_tab()
    {

        int m_l = LIGNES/2+1;
        int m_c = COL/2+1;
        for (size_t u(0); u<LIGNES; u++)
        {
            for(size_t v(0); v<COL; v++)
            {
                display_tab[u][v] = ' ';
                display_tab_color[u][v] = 15;
            }
        }
        for(size_t i(0); i<list_objects.size(); i++)
        {
            if (list_objects[i]->m_name != "Camera")//ne pas afficher la caméra ^^
            {
                size_t nb_points = list_objects[i]->m_definig_points.size();
                vector<Point> list_points = list_objects[i]->m_definig_points;
                double length = list_objects[i]->m_length;
                int color = list_objects_color[i];

                for(size_t j(0); j<nb_points; j++)
                {
                    double x(0);
                    double y(0);
                    double z(0);
                    if (is_camera_view)
                    {
                        x = list_points[j].m_x + list_objects_coord[i].m_x;
                        y = list_points[j].m_y + list_objects_coord[i].m_y;
                        z = list_points[j].m_z + list_objects_coord[i].m_z;
                        Point temp(m_cam.transformation_plan(Point(x,y,z)));
                        x = temp.m_x;
                        y = temp.m_y;
                        z = temp.m_z;
                        //cout<<"JE PAAAASSSSEEE"<<Sendl;
                    }
                    else
                    {
                        x = list_points[j].m_x + list_objects_coord[i].m_x;
                        y = list_points[j].m_y + list_objects_coord[i].m_y;
                        z = list_points[j].m_z + list_objects_coord[i].m_z;
                    }

                    int r_y = round(y);
                    int r_z = round(z);
                    int r_z_axis = round(m_z_axis);
                    int r_y_axis = round(m_y_axis);
                    //cout<<"x: "<<r_x<<endl;
                    //cout<<"y: "<<r_y<<endl;
                    //cout<<"z: "<<r_z<<endl;
                    //cout<<"m_c : "<<-m_c<<endl;
                    //cout<<"m_l : "<<-m_l<<endl;


                    ///revoir le concept
                    //Camera cam_copy = get_camera();
                    //if(r_z < m_c && r_z > -m_c && r_y < m_c && r_y > -m_l)//dans le tableau
                    if(m_l+r_z+r_z_axis < LIGNES && m_c+r_y+r_y_axis < COL)
                    {
                        //cout<<"dans le tableau"<<endl;
                        if(x >= -sqrt(pow(length,2)+pow(length,2)) && x < -sqrt(pow(length,2)+pow(length,2))/2)
                        {
                            if(display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] != '*' && display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] != '#' && display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] != '@')

                            {
                                display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = '.';
                                display_tab_color[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = color;
                            }
                        }
                        else if(x >= -sqrt(pow(length,2)+pow(length,2))/2 && x < 0)
                        {
                            if(display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] != '#' && display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] != '@')
                            {
                                display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = '*';
                                display_tab_color[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = color;
                            }

                        }
                        else if(x >= 0 && x < sqrt(pow(length,2)+pow(length,2))/2)
                        {
                            if(display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] != '@')
                            {
                                display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = '#';
                                display_tab_color[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = color;
                            }
                        }
                        else if(x >= sqrt(pow(length,2)+pow(length,2))/2 && x <= sqrt(pow(length,2)+pow(length,2)))
                        {
                            display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = '@';
                            display_tab_color[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = color;
                        }
                        else
                        {
                            display_tab[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = 49+j;
                            display_tab_color[LIGNES-m_l-r_z-r_z_axis][COL-m_c-r_y-r_y_axis] = color;
                        }
                   }

                }
            }
        }
    }
    void display()
    {
        string tableau_fini;
        system("CLS");
        set_display_tab();
        if(is_camera_view)
        {
            set_point_color(15);
            cout<<"Camera"<<endl;
        }
        else
        {
            set_point_color(list_objects_color[m_current_object_index]);
            cout<<list_objects[m_current_object_index]->m_name<<endl;
        }

        reset_text_color();
        for(size_t i(0); i<LIGNES; i++)
        {
            for(size_t j(0); j<COL; j++)
            {
                set_point_color(display_tab_color[i][j]);
                cout<<display_tab[i][j];
                //tableau_fini+=display_tab[i][j];
                //cout<<;
                //reset_text_color();
            }
            cout<<i<<"\n"<<endl;
            //tableau_fini+=(to_string(i)+"\n");
        }
        //cout<<tableau_fini;
    }
};

int main()
{
    Environnement E(0,0,0);

    Object sphere;
    sphere.make_sphere(10,0.2,"sphere");

    Object cube;
    cube.make_cube(10, 0.2, "cube");

    Object line;
    line.make_line(12,1,"ligne");
    //E.add_object(Camera,Point(0,0,0), "blanc")
    E.add_object(&sphere, Point(0,0,0), "rouge");
    E.add_object(&cube, Point(0,0,0), "bleu");
    E.add_object(&line, Point(0,0,0), "jaune");
    E.display();
    while(true)
    {
        if(key_pressed(E))
        {
            E.display();
        }
    }

    return 0;
}

Arrete point_to_arrrete(Point p)
{
    return Arrete(p.m_x, p.m_y, p.m_z);
}

Point arrete_to_point(Arrete a)
{
    return Point(a.m_x, a.m_y, a.m_z);
}

bool key_pressed(Environnement &E)
{
    bool is_a_key_pressed(true);

    if(GetKeyState(VK_TAB) & 0x8000)
    {
        E.select_next_object();
    }

    else if(GetKeyState('C') & 0x8000)
    {
        if (E.is_camera_view)
        {
            E.is_camera_view = false;
        }
        else
        {
            E.is_camera_view = true;
        }
    }

    else if(GetKeyState('Q') & 0x8000)
    {
        if (E.is_camera_view)
        {
            E.m_cam.rotate_z(0.1);
        }
        else
        {
            E.list_objects[E.m_current_object_index]->rotate_x(0.1);
        }
    }
    else if(GetKeyState('E') & 0x8000)
    {
        if (E.is_camera_view)
        {
            E.m_cam.rotate_z(-0.1);
        }
        else
        {
            E.list_objects[E.m_current_object_index]->rotate_x(-0.1);
        }
    }
    else if(GetKeyState('W') & 0x8000)
    {
        if (E.is_camera_view)
        {
            E.m_cam.translate({1,0,0});
        }
        else
        {
            E.list_objects[E.m_current_object_index]->rotate_y(0.1);
        }

    }
    else if(GetKeyState('S') & 0x8000)
    {
        if (E.is_camera_view)
        {
            E.m_cam.translate({1,0,0});
        }
        else
        {
            E.list_objects[E.m_current_object_index]->rotate_y(-0.1);
        }
    }
    else if(GetKeyState('A') & 0x8000)//is down
    {
        if (E.is_camera_view)
        {
            E.m_cam.translate({0,-1,0});
        }
        else
        {
            E.list_objects[E.m_current_object_index]->rotate_z(-0.1);
        }
    }
    else if(GetKeyState('D') & 0x8000)
    {
        if (E.is_camera_view)
        {
            E.m_cam.translate({0,1,0});
        }
        else
        {
            E.list_objects[E.m_current_object_index]->rotate_z(0.1);
        }
    }
    else if(GetKeyState(VK_UP) & 0x8000)
    {
        Point p = Point(0,0,-1);
        E.move_object(p);
    }
    else if(GetKeyState(VK_DOWN) & 0x8000)
    {
        Point p = Point(0,0,1);
        E.move_object(p);
    }
    else if(GetKeyState(VK_LEFT) & 0x8000)
    {
        Point p = Point(0,-1,0);
        E.move_object(p);
    }
    else if(GetKeyState(VK_RIGHT) & 0x8000)
    {
        Point p = Point(0,1,0);
        E.move_object(p);
    }
    else if(GetKeyState('I') & 0x8000)
    {
        E.m_z_axis -= 1;
    }
    else if(GetKeyState('K') & 0x8000)
    {
        E.m_z_axis += 1;
    }
    else if(GetKeyState('J') & 0x8000)
    {
        E.m_y_axis -= 1;
    }
    else if(GetKeyState('L') & 0x8000)
    {
        E.m_y_axis += 1;
    }
    else if(GetKeyState(VK_ESCAPE) & 0x8000)
    {
        afficher_texte_color("***MERCI DE VOTRE VISITE***", "vert");
        getch();

        exit(0);
    }
    else
        is_a_key_pressed = false;

    return is_a_key_pressed;
}


