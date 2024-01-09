#define _USE_MATH_DEFINES
#include <chrono>
#include "M_matrices.h"
#include "Color_picker.h"
#include <Windows.h>
#include <conio.h>

///probleme avec les axes : axes inversés
using namespace std;

struct Point;
struct Arrete;
struct Object;
struct Environnement;
Arrete point_to_arrrete(Point p);
Point arrete_to_point(Arrete a);
bool key_pressed(Environnement &E);

const size_t LIGNES(40), COL(50);//taille de l'environnement

struct Point
{
    Point(){}
    Point(double x, double y, double z): m_x(x), m_y(y), m_z(z){}
    double m_x;
    double m_y;
    double m_z;
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
};


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
Point operator*(Point a, double b)
{
    Point S;
    S.m_x = a.m_x * b;
    S.m_y = a.m_y * b;
    S.m_z = a.m_z * b;
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
    double m_rot_x;
    double m_rot_y;
    double m_rot_z;
    double m_length;
    double m_rayon;
    double m_pas;
    M_matrices m_matrice_rot = get_matrice_identite(3,3);
    Point m_origine;
    vector<Point> m_definig_points;
    vector<Arrete> m_arretes;
    string m_name;

    void make_line(double length, double pas, string name)
    {
        m_name = name;
        m_definig_points.clear();
        m_length = length;
        m_pas = pas;

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

        for(int i = 0; i*pas<rayon; i++)//de 0 � 2*rayon
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
        //derni�res arr�tes
        a = point_to_arrrete(m_definig_points[3]);//1 � 4
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[0];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[5]);//1 � 6
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[0];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[6] - m_definig_points[1]);//2 � 7
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[1];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[7] - m_definig_points[2]);//3 � 8
        m_arretes.push_back(a);
        for(size_t j(1) ; j*pas<length; j++)
            {
                Point o = m_definig_points[2];
                a.transform_to_unit_vect();
                Arrete v = a *(j*pas);///peut manquer une arrete
                m_definig_points.push_back(o+arrete_to_point(v));
            }
        a = point_to_arrrete(m_definig_points[7] - m_definig_points[4]);//5 � 8
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

        set_centered();
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
            for(size_t i(0); i<m_definig_points.size(); i++)
            {
                Point temp = Point(translation[0], translation[1], translation[2]);
                translate_point(m_definig_points[i], temp);
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

struct Environnement
{
    Environnement():m_x_axis(0), m_y_axis(0), m_z_axis(0), m_x_vect(Point(1,0,0)), m_y_vect(Point(0,1,0)), m_z_vect(Point(0,0,1))
    {
    }
    Environnement(double x_axis, double y_axis, double z_axis): m_x_axis(x_axis), m_y_axis(y_axis), m_z_axis(z_axis), m_x_vect(Point(1,0,0)), m_y_vect(Point(0,1,0)), m_z_vect(Point(0,0,1)){}
    //static const size_t LIGNES=40, COL=50;
    char display_tab[LIGNES][COL];
    int display_tab_color[LIGNES][COL];
    int m_current_object_index = 0;

    M_matrices m_matrice_rot = get_matrice_identite(3,3);
    double m_x_axis;//coord des axes
    double m_y_axis;
    double m_z_axis;
    double m_length_axis;
    double m_pas_axis;
    Point m_x_vect, m_y_vect, m_z_vect;
    vector<Point> axis_vect;
    vector<Point> axis_vect_definig_points;
    vector<Object> list_objects;
    vector<int> list_objects_color;
    vector<Point> list_objects_coord;
    //vector<string> COLor_pannel = {"bleu","rouge"};
    void rotate_axis()
    {
        M_matrices x_coord(3,1), y_coord(3,1), z_coord(3,1);
        x_coord.set_matrice_vect(m_x_vect.point_to_vector());
        y_coord.set_matrice_vect(m_y_vect.point_to_vector());
        z_coord.set_matrice_vect(m_z_vect.point_to_vector());
        M_matrices temp_x = m_matrice_rot*x_coord;
        M_matrices temp_y = m_matrice_rot*y_coord;
        M_matrices temp_z = m_matrice_rot*z_coord;

        m_x_vect.set_coord(temp_x.get_matrice_vect());
        m_y_vect.set_coord(temp_y.get_matrice_vect());
        m_z_vect.set_coord(temp_z.get_matrice_vect());
    }
    void rotate_x(double angle_rad)
    {
        vector<double> matrice_rot_x = {1,0,0,
                                        0,cos(angle_rad),-sin(angle_rad),
                                        0,sin(angle_rad),cos(angle_rad)};

        m_matrice_rot.set_matrice_vect(matrice_rot_x);
        rotate_axis();

    }
    void rotate_y(double angle_rad)
    {
        vector<double> matrice_rot_y = {cos(angle_rad),0,sin(angle_rad),
                                        0,1,0,
                                        -sin(angle_rad),0,cos(angle_rad)};
        m_matrice_rot.set_matrice_vect(matrice_rot_y);
        rotate_axis();
    }
    void rotate_z(double angle_rad)
    {
        vector<double> matrice_rot_z = {cos(angle_rad),-sin(angle_rad),0,
                                        sin(angle_rad),cos(angle_rad) ,0,
                                        0             ,0              ,1};
        m_matrice_rot.set_matrice_vect(matrice_rot_z);
        rotate_axis();
    }

    void select_next_object()
    {
        if (m_current_object_index == list_objects.size()-1)
            m_current_object_index = 0;
        else
            m_current_object_index++;
    }

    void move_object(Point target)//d�place l'objet courant dans l'environnement
    {
        list_objects_coord[m_current_object_index] += target;
    }

    void clear_all()
    {
        list_objects.clear();
        list_objects_color.clear();
        list_objects_coord.clear();
    }
    void add_object(Object item, Point coord, string color)
    {
        list_objects.push_back(item);
        int color_int = transform_text_color(color);
        list_objects_color.push_back(color_int);
        list_objects_coord.push_back(coord);
    }
    void make_axis(double pas, double length)
    {
        m_matrice_rot.set_nb_lignes_col(3,3);
        m_length_axis = length;
        m_pas_axis = pas;
        axis_vect.push_back(m_x_vect);
        axis_vect.push_back(m_y_vect);
        axis_vect.push_back(m_z_vect);
        Point o = Point(m_x_axis, m_y_axis, m_z_axis);
        axis_vect_definig_points.push_back(o);//Point d'intersection des 3 axes
        for(size_t i(0); i<axis_vect.size(); i++)
        {
            for(size_t j(1); j*pas<length; j++)
            {
                Point p;
                p = axis_vect[i]*(j*pas);
                axis_vect_definig_points.push_back(o+p);
            }
        }
    }
    void refresh_axis()
    {
        axis_vect_definig_points.clear();
        axis_vect[0] = m_x_vect;
        axis_vect[1] = m_y_vect;
        axis_vect[2] = m_z_vect;
        Point o = Point(m_x_axis, m_y_axis, m_z_axis);
        axis_vect_definig_points.push_back(o);//Point d'intersection des 3 axes
        for(size_t i(0); i<axis_vect.size(); i++)
        {
            for(size_t j(1); j*m_pas_axis<m_length_axis; j++)
            {
                Point p;
                p = axis_vect[i]*(j*m_pas_axis);
                axis_vect_definig_points.push_back(o+p);
                //cout<<p.m_x<<" "<<p.m_y<<" "<<p.m_z<<endl;
            }
        }
        cout<<m_x_vect.m_x<<" "<<m_x_vect.m_y<<" "<<m_x_vect.m_z<<endl;
    }

    void set_environnement_coord(double x_axis, double y_axis, double z_axis)
    {
        m_x_axis = x_axis;
        m_y_axis = y_axis;
        m_z_axis = z_axis;
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

        size_t m_l = LIGNES/2;
        size_t m_c = COL/2;
        clear_display_tab();
        for(size_t i(0); i<list_objects.size(); i++)
        {
            size_t nb_points = list_objects[i].m_definig_points.size();
            vector<Point> list_points = list_objects[i].m_definig_points;
            double length = list_objects[i].m_length;
            int color = list_objects_color[i];

            for(size_t j(0); j<nb_points; j++)
            {
                double x = list_points[j].m_x + list_objects_coord[i].m_x + m_x_axis;
                double y = list_points[j].m_y + list_objects_coord[i].m_y + m_y_axis;
                double z = list_points[j].m_z + list_objects_coord[i].m_z + m_z_axis;
                int r_y = round(y);
                int r_z = round(z);
                //int r_z_axis = round(m_z_axis);
                //int r_y_axis = round(m_y_axis);
                if(list_points[j].m_x >= -sqrt(pow(length,2)+pow(length,2)) && list_points[j].m_x < -sqrt(pow(length,2)+pow(length,2))/2)//marche pour la taille du cube
                {
                    if(display_tab[m_c+r_y][m_l-r_z] != '*' && display_tab[25+r_y][25-r_z] != '#' && display_tab[25+r_y][25-r_z] != '@')
                    {
                        display_tab[m_c+r_y][m_l-r_z] = '.';
                        display_tab_color[m_c+r_y][m_l-r_z] = color;
                    }
                }
                else if(list_points[j].m_x >= -sqrt(pow(length,2)+pow(length,2))/2 && list_points[j].m_x < 0)
                {
                    if(display_tab[m_c+r_y][m_l-r_z] != '#' && display_tab[25-r_z][25+r_y] != '@')
                    {
                        display_tab[m_c+r_y][m_l-r_z] = '*';
                        display_tab_color[m_c+r_y][m_l-r_z] = color;
                    }

                }
                else if(list_points[j].m_x >= 0 && list_points[j].m_x < sqrt(pow(length,2)+pow(length,2))/2)
                {
                    if(display_tab[m_c+r_y][m_l-r_z] != '@')
                    {
                        display_tab[m_c+r_y][m_l-r_z] = '#';
                        display_tab_color[m_c+r_y][m_l-r_z] = color;
                    }
                }
                else if(list_points[j].m_x >= sqrt(pow(length,2)+pow(length,2))/2 && list_points[j].m_x <= sqrt(pow(length,2)+pow(length,2)))
                {
                    display_tab[m_c+r_y][m_l-r_z] = '@';
                    display_tab_color[m_c+r_y][m_l-r_z] = color;
                }
                else
                {
                    display_tab[m_c+r_y][m_l-r_z] = 49+j;
                    display_tab_color[m_c+r_y][m_l-r_z] = color;
                }
            }
        }

        Point p = axis_vect_definig_points[0];
        int r_axis_y = round(p.m_y), r_axis_z = round(p.m_z);
        display_tab[r_axis_y][r_axis_z] = 'O';
        display_tab_color[r_axis_y][r_axis_z] = 13;//violet
        for(size_t i(1); i<axis_vect_definig_points.size(); i++)/////////////Probleme
        {
            //cout<<"mdr"<<i<<endl;
            p = axis_vect_definig_points[i];
            r_axis_y = round(p.m_y);
            r_axis_z = round(p.m_z);
            if(i<=(axis_vect_definig_points.size()-1)/3)
            {
                display_tab_color[r_axis_y][r_axis_z] = 4;//rouge x
                if(i == (axis_vect_definig_points.size()-1)/3)
                    display_tab[r_axis_y][r_axis_z] = 'X';
                else
                    display_tab[r_axis_y][r_axis_z] = '%';
            }
            else if(i<=2*(axis_vect_definig_points.size()-1)/3)
            {
                display_tab_color[r_axis_y][r_axis_z] = 10;//vert y
                if(i == 2*(axis_vect_definig_points.size()-1)/3)
                    display_tab[r_axis_y][r_axis_z] = 'Y';
                   // cout<<"ok"<<endl;
                else
                    display_tab[r_axis_y][r_axis_z] = '%';
            }
            else
            {
                display_tab_color[r_axis_y][r_axis_z] = 9;//bleu z
                if(i == (axis_vect_definig_points.size()-1))
                    display_tab[r_axis_y][r_axis_z] = 'Z';
                else
                    display_tab[r_axis_y][r_axis_z] = '%';
            }
        }
    }
    void display()
    {
        system("CLS");
        set_display_tab();
        set_point_color(list_objects_color[m_current_object_index]);
        cout<<list_objects[m_current_object_index].m_name<<endl;
        reset_text_color();
        for(size_t i(0); i<LIGNES; i++)
        {
            for(size_t j(0); j<COL; j++)
            {
                set_point_color(display_tab_color[i][j]);
                cout<<display_tab[i][j];
                reset_text_color();
            }
            cout<<i<<endl;
        }

    }
};

int main()
{
    Environnement E;
    E.make_axis(1,15);

    Object sphere;
    sphere.make_sphere(10,0.2,"sphere");

    Object cube;
    cube.make_cube(12, 0.2, "cube");

    Object line;
    line.make_line(12,1,"ligne");

    E.add_object(sphere, Point(0,0,0), "rouge");
    E.add_object(cube, Point(0,0,0), "bleu");
    E.add_object(line, Point(0,0,0), "jaune");
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

    else if(GetKeyState('Q') & 0x8000)
    {
        E.list_objects[E.m_current_object_index].rotate_x(0.1);
    }
    else if(GetKeyState('E') & 0x8000)
    {
        E.list_objects[E.m_current_object_index].rotate_x(-0.1);
    }
    else if(GetKeyState('W') & 0x8000)
    {
        E.list_objects[E.m_current_object_index].rotate_z(-0.1);
    }
    else if(GetKeyState('S') & 0x8000)
    {
        E.list_objects[E.m_current_object_index].rotate_z(0.1);
    }
    else if(GetKeyState('A') & 0x8000)//is down
    {
        E.list_objects[E.m_current_object_index].rotate_y(-0.1);
    }
    else if(GetKeyState('D') & 0x8000)
    {
        E.list_objects[E.m_current_object_index].rotate_y(0.1);
    }
    else if(GetKeyState(VK_UP) & 0x8000)
    {
        Point p = Point(0,-1,0);
        E.move_object(p);
    }
    else if(GetKeyState(VK_DOWN) & 0x8000)
    {
        Point p = Point(0,1,0);
        E.move_object(p);
    }
    else if(GetKeyState(VK_LEFT) & 0x8000)
    {
        Point p = Point(0,0,1);
        E.move_object(p);
    }
    else if(GetKeyState(VK_RIGHT) & 0x8000)
    {
        Point p = Point(0,0,-1);
        E.move_object(p);
    }
    else if(GetKeyState('I') & 0x8000)
    {
        E.m_y_axis -= 1;
        E.refresh_axis();
    }
    else if(GetKeyState('K') & 0x8000)
    {
        E.m_y_axis += 1;
        E.refresh_axis();
    }
    else if(GetKeyState('J') & 0x8000)
    {
        E.m_z_axis -= 1;
        E.refresh_axis();
    }
    else if(GetKeyState('L') & 0x8000)
    {
        E.m_z_axis += 1;
        E.refresh_axis();
    }
    else if(GetKeyState('1') & 0x8000)
    {
        E.rotate_x(0.2);
        E.refresh_axis();
    }
    else if(GetKeyState('2') & 0x8000)
    {
        E.rotate_y(0.2);
        E.refresh_axis();
    }
    else if(GetKeyState('3') & 0x8000)
    {
        E.rotate_z(0.2);
        E.refresh_axis();
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


