/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
 public class Edge
 {
    public int start;
    public int end;
    
    public Edge()
    {
        this.start =-1;
        this.end = -1;
    }
    
    public Edge(int s, int e)
    {
        this.start = s;
        this.end = e;
    }
    
    public Edge(Edge e)
    {
        this.start = e.start;
        this.end = e.end;
    }
}
