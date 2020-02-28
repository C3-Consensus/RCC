/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */

public class Atom{
    public int id;
    public String name;
    public Point3d coords;
    
    public Atom(int id, String name)
    {
        this.id = id;
        this.name = new String(name);
        coords = new Point3d();
    }

    public Atom(int id, String name, Point3d a) 
    {
        this.id = id;
        this.name = new String(name);
        coords = new Point3d(a);
    }
    public Atom(int id, String name, double x, double y, double z) 
    {
        this.id = id;
        this.name = new String(name);
        this.coords = new Point3d(x,y,z);
    }
    
    public Atom(Atom a)
    {
        this.id = a.id;
        this.name = new String(a.name);
        this.coords = new Point3d(a.coords);
    }
    
}
