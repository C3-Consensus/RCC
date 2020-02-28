/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Backbone {
    public Atom N;
    public Atom Ca;
    public Atom C;
    public Atom O;
    
    public Backbone()
    {
        N = null;
        Ca = null;
        C = null;
        O = null;
    }
    
    public Backbone(Atom N, Atom Ca, Atom C, Atom O)
    {
        this.N = N;
        this.C = C;
        this.Ca = Ca;
        this.O = O;
    }
    
    public Backbone(Backbone a)
    {
        this.N = null;
        this.C = null;
        this.Ca = null;
        this.O = null;
        
        if(a.N != null)
        {
            this.N = new Atom(a.N);
        }
        
        if(a.Ca != null)
        {
            this.Ca = new Atom(a.Ca);            
        }
        if(a.C != null)
        {
            this.C = new Atom(a.C);
        }
        if(a.O != null)
        {
            this.O = new Atom(a.O);
        }
    }
}
