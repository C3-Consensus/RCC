
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class CATHDomain {
    
    static private String ProcessDomain(Peptide fullPeptide, Domain dom, int numDomain, String fileName, int totalDomains) throws IOException, FileNotFoundException, IllegalArgumentException, IllegalAccessException {
        for(int i=0; i< dom.numSegments; i++)
        {
            if( !dom.chainStart[0].equalsIgnoreCase( dom.chainStart[i]) || !dom.chainStart[0].equalsIgnoreCase( dom.chainEnd[i]) )
            {
                System.err.println("Different chain, file "+ fileName + " domain " + numDomain );
            }
        }
        
        if(! fullPeptide.chain.equalsIgnoreCase( dom.chainStart[0] ) )
        {
            System.err.println("Chain from domain does not correspond to file");
        }
        
        Peptide p = fullPeptide.ExtractDomain(dom);
        String graphFileName = Parameters.pathGraphs+fileName+p.chain+"D"+numDomain+".csv";
        if(numDomain < 10)
        {
            dom.name = fileName+p.chain+"0"+ (numDomain+1);
        }
        else
        {
            dom.name = fileName+p.chain + (numDomain+1);
        }
        
        if(totalDomains == 1)
        {
            dom.name = fileName+p.chain+"00";
        }
        
        p.buildGraph(Parameters.min, Parameters.max, Parameters.distance_criterium, Parameters.pairs, Parameters.flagBackbone, graphFileName);
        Rcc r = new Rcc(graphFileName,"tomita",null);
        String info = r.printable;
        //String info = GetRCCFromPDB.printRCC(graphFileName, "tomita");
        
        info += "\t"+p.residues.size() + "\t" + p.missingResidues;
        
        if( p.missingResidues != 0)
        {
            info += "\t"+p.comments;
        }
        
        return info;
    }
    
    static void ProcessFile(String file, String outputFileName, String errorFileName) throws FileNotFoundException, IOException, IllegalArgumentException, IllegalAccessException
    {
        BufferedReader inFile = new BufferedReader(new FileReader(file));
        BufferedWriter outFile = new BufferedWriter(new FileWriter(outputFileName));
        BufferedWriter errorFile = new BufferedWriter(new FileWriter(errorFileName));
        

        String line;
        int linesProcessed = 0;
        while ( (line = inFile.readLine())!= null)
        {
            if(line.startsWith("#"))
            {
                continue;
            }
                        
            linesProcessed ++;
            if( linesProcessed%1000 == 0)
            {
                System.out.println("lines processed:" + linesProcessed );
                System.out.println("current line:" + line );
            }

/*
0      1    2   3  4    5 6 7   8  9  10 11 12 13 14 15 16 17  18 19 20 21  22
1cnsA  D02 F00  2  A    1 - A   87 -  A  146 - A  243 -  1  A   88 - A  145 -
                N |C    S I C    E I| C    S I C    E I| N |C    S I C    E I|
               |<--------------Domain One------------->|<-----Domain Two---->|
                  |<--Segment One-->|<---Segment Two-->|   |<--Segment One-->|            
*/            
            //parse line
            String[] fields = line.split("\\s+");
            String fileName = fields[0].substring(0, fields[0].length()-1);
            String currentChain = "" + fields[0].charAt( fields[0].length()-1 );
            String subFolder = fileName.substring(1,3);
            
            //load protein
            String pdbFileName = Parameters.pathPdb+ subFolder +"/pdb"+  fileName+".ent.gz";
            
            //skip files missing
            File f = new File(pdbFileName);
            if( !f.exists()) 
            {
                System.out.println("File not found: "+ pdbFileName);
                errorFile.write("file not fount:"+pdbFileName+"\n");
                continue;
            } 
            
            Peptide p = new Peptide( pdbFileName, currentChain, Parameters.flagBackbone);
            
            int numDomains = Integer.parseInt(fields[1].substring(1));
            int nextField = 3;
            for(int domain = 0; domain< numDomains; domain++)
            {
                int numSegments = Integer.parseInt(fields[nextField++]);
                Domain dom = new Domain(numSegments);
                //extract segments
                for(int segment=0; segment < numSegments; segment++)
                {
                    // Chain Start Insert Chain End Insert
                    dom.chainStart[segment] = fields[nextField++];
                    dom.start[segment] = Integer.parseInt(fields[nextField++]);
                    if( !fields[nextField++].equalsIgnoreCase("-"))
                    {
                        System.err.println( "Has insert character\n" + line);
                    }
                    
                    dom.chainEnd[segment] = fields[nextField++];
                    dom.end[segment] = Integer.parseInt(fields[nextField++]);
                    if( !fields[nextField++].equalsIgnoreCase("-"))
                    {
                        System.err.println( "Has insert character\n" + line);
                    }                    
                }                                
                String rccMessage = ProcessDomain(p, dom,domain,fileName, numDomains);
                outFile.write(dom.name + "\t" + rccMessage+"\n");
            }
        }
        inFile.close();
        outFile.close();
        errorFile.close();
    }
}
