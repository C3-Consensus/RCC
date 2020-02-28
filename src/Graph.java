
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Graph {
    ArrayList<Edge> edges;
    public Graph ()
    {
        edges = new ArrayList<>();
    }

    public void printGraph (String filename, int numResidues) throws IOException
    {
        StringBuffer bufferS = new StringBuffer();
        BufferedWriter outfile = new BufferedWriter(new FileWriter(filename));
        int order = 0;
        for(int e=0; e < edges.size(); e++)
        {
            bufferS.append( edges.get(e).start + "," + edges.get(e).end);
            bufferS.append( "\n" + edges.get(e).end + "," + edges.get(e).start+"\n");
            order += 1;
        }

        outfile.write( numResidues +"\n"+(2*order)+"\n");
        outfile.write(bufferS.toString());
        outfile.flush();
        outfile.close();

    }
}
