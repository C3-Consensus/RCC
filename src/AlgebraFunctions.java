/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class AlgebraFunctions {
    public static Point3d e = new Point3d(1, 1, 1);
    
    public static double calculateDistance(double x1, double y1, double z1, double x2, double y2, double z2)
    {
        return Math.sqrt(  Math.pow(x1-x2, 2)   +   Math.pow(y1-y2, 2)  +   Math.pow(z1-z2, 2)  );
    }

    public static double calculateDistance(Point3d a, Point3d b)
    {
        return calculateDistance(a.x, a.y, a.z, b.x, b.y, b.z);        
    }
    public static double volumeIntersection(double d, double r1, double r2)
    {
        //http://mathworld.wolfram.com/Sphere-SphereIntersection.html
        // V = 	(pi(R+r-d)^2(d^2+2dr-3r^2+2dR+6rR-3R^2))/(12d).
        return Math.PI * 
                    ( Math.pow( r1 + r2 - d , 2 ) 
                        * ( Math.pow(d, 2) 
                            + 2 * d * (r1 + r2) 
                            + 6*r1*r2 
                            - 3*Math.pow(r1, 2) 
                            - 3*Math.pow(r2, 2) 
                        )
                    )
                / (12.0 * d);
    }	
    
    public static double CalculateAngle(Point3d a, Point3d b, Point3d axis)
    {
        //     dot = x1*x2        + y1*y2        + z1*z2
        double dot = (a.x * b.x ) + (a.y * b.y ) + (a.z * b.z);
        //     det = x1*y2*zn       + x2*yn*z1       + xn*y1*z2       - z1*y2*xn       - z2*yn*x1       - zn*y1*x2
        double det = a.x*b.y*axis.z + b.x*axis.y*a.z + axis.x*a.y*b.z - a.z*b.y*axis.x - b.z*axis.y*a.x - axis.z*a.y*b.x;
        //angle = atan2(det, dot)
        return Math.atan2(det,dot);
    }
    
    public static double InnerProduct( Point3d a, Point3d b)
    {
        return (a.x * b.x ) + (a.y * b.y ) + (a.z * b.z);
    }
    
    public static Point3d SubstractPoints( Point3d a, Point3d b)
    {
        Point3d c = new Point3d();
        c.x = a.x - b.x;
        c.y = a.y - b.y;
        c.z = a.z - b.z;
        return c;
    }
    
    public static Point3d CalculateNormalVector( Point3d u, Point3d v)
    {
        Point3d n = new Point3d();
        
        n.x = e.x - u.x - v.x;
        n.y = e.y - u.y - v.y;
        n.z = e.z - u.z - v.z;
        return n;
    }

    public static Point3d CrossProduct(Point3d u, Point3d v)
    {
        return new Point3d(
    //s_x = u_y*v_z - u_z*v_y
            u.y*v.z - u.z*v.y,
    //s_y = u_z*v_x - u_x*v_z
            u.z*v.x - u.x*v.z,
    //s_z = u_x*v_y - u_y*v_x
            u.x*v.y - u.y*v.x
        );
    }
    public static Point3d CrossProductNormalized(Point3d u, Point3d v)
    {
        Point3d n = CrossProduct(u, v);
        n.normalize();
        return n;
    }
}
