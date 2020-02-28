/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Point3d {
    public double x;
    public double y;
    public double z;    
        
    public Point3d(double x, double y, double z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    Point3d() 
    {
    }
    
    public Point3d(Point3d p, Point3d q)
    {
        this.x = p.x - q.x;
        this.y = p.y - q.y;
        this.z = p.z - q.z;
    }
    
    public Point3d(Point3d p)
    {
        this.x = p.x;
        this.y = p.y;
        this.z = p.z;
    }
    
    public void setCoords(double x, double y, double z)
    {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public void multiplyByScalar(double lamda)
    {
        this.x *= lamda;
        this.y *= lamda;
        this.z *= lamda;        
    }
    
    public Point3d getInverse()
    {
        return new Point3d(-this.x, -this.y, -this.z);
    }
    
    public double getNorm()
    {
        return Math.sqrt( (this.x*this.x) + (this.y * this.y) + (this.z*this.z)  );
    }
    
    public Point3d getScaled( double lamda)
    {
        return new Point3d( this.x* lamda, this.y*lamda, this.z*lamda);
    }
    
    public void normalize()
    {
        double n = 1.0 / this.getNorm();
        this.x *= n;
        this.y *= n;
        this.z *= n;
    }
    
    public void sumVector(Point3d v)
    {
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;
    }
    
    public void sumAtom(Atom v)
    {
        this.x += v.coords.x;
        this.y += v.coords.y;
        this.z += v.coords.z;
    }
    
    public void rotateWithRodriguez(double theta, Point3d axis)
    {
    // v_r = v cos (t) +(k x v) sin(t) + k < k, v> (1 - cos (t) )
        Point3d KxV = AlgebraFunctions.CrossProduct(this, axis);
        double KV = AlgebraFunctions.InnerProduct(this, axis);
        double sint = Math.sin(theta);
        double cost = Math.cos(theta);
        this.x = (this.x * cost) + (KxV.x * sint) + ( axis.x * KV * (1-cost) );
        this.y = (this.y * cost) + (KxV.y * sint) + ( axis.y * KV * (1-cost) );
        this.z = (this.z * cost) + (KxV.z * sint) + ( axis.z * KV * (1-cost) );        
    }
}
