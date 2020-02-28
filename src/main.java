
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Time;
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



public class main {
    
    public static void main(String[] args) throws IOException, FileNotFoundException, IllegalArgumentException, IllegalAccessException, InterruptedException 
    {

        Parameters.LoadParameters(args);
        
        if(!Parameters.loadCorrect)
        {
            return;
        }
        
        if(Parameters.fileList != null)
        {
            BufferedReader fileList = new BufferedReader(new FileReader(Parameters.fileList));
            String line = "";
            
            while( (line = fileList.readLine())!= null)
            {
                String fields[] = line.split(",");
                Parameters.pdbFileName = fields[0];
                Parameters.chain = fields[1];
                
                String localFilePath = Parameters.pathPdb + Parameters.pdbFileName.substring(1, 3).toLowerCase()+"/pdb"+Parameters.pdbFileName.toLowerCase()+".ent.gz";
                //String localFilePath = Parameters.pathPdb +Parameters.pdbFileName.toLowerCase()+".cif.gz";
                
                try
                {
                    Peptide p = new Peptide( localFilePath, Parameters.chain, Parameters.flagLateralChain);
                    if(p.residues.size() == 0)
                    {
                        System.out.println("file " + Parameters.pdbFileName + "  chain " + Parameters.chain + " is empty");
                        continue;
                    }

                    Parameters.graphFileName = fields[0]+"-graph.txt";
                    p.buildGraph(Parameters.min, Parameters.max, Parameters.distance_criterium, Parameters.pairs, Parameters.flagLateralChain, Parameters.graphFileName);
                    Rcc r = new Rcc(Parameters.graphFileName,"tomita",Parameters.pdbFileName+"cliques.txt");
                    System.out.println(line+","+r.printable);
                }
                catch (java.io.FileNotFoundException e)
                {
                    System.out.println("Not found: " + fields[0] + "\tat\t"+localFilePath);
                }
                
            }
            return;
        }
        
        
        Peptide p = new Peptide (Parameters.pdbFileName, Parameters.chain, Parameters.flagLateralChain);

        if(p.residues.size() == 0)
        {
            System.out.println("file " + Parameters.pdbFileName + "  chain " + Parameters.chain + " is empty");
        }

        p.buildGraph(Parameters.min, Parameters.max, Parameters.distance_criterium, Parameters.pairs, Parameters.flagLateralChain, Parameters.graphFileName);
        Rcc r = new Rcc(Parameters.graphFileName,"tomita",Parameters.pdbFileName+"cliques.txt");
        
        System.out.println(r.printable);
        
        return;
        
    }

}
