
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class Clique {
    public int[] nodes;
    public int start;
    public int end;
    //public String typeS;
    public int rccEntry;
    
    public Clique(Clique c)
    {
        this.start = c.start;
        this.end = c.end;
        this.rccEntry = rccEntry;
        this.nodes = new int[c.nodes.length];
        for(int n = 0; n< c.nodes.length; n++)
        {
            this.nodes[n] = c.nodes[n];
        }
    }

    public Clique(ArrayList<Integer> nodes, int rcc)
    {
        this.nodes = new int[nodes.size()];
        start = (int) Constants.infPos;
        end = (int) Constants.infNeg;
        int j=0;
        for(Integer i : nodes)
        {
            this.nodes[j] = i;
            j++;
            start = Math.min(start, i);
            end = Math.max(end, i);
        }
        this.rccEntry = rcc;
    }
}
