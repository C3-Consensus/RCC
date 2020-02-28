
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class Rcc {
    ArrayList<Integer> rcc;
    String printable;
    ArrayList<Clique> cliques;
    
    public Rcc(Rcc r)
    {
        this.rcc = new  ArrayList<>();
        for(Integer i: r.rcc)
        {
            this.rcc.add(new Integer(i));
        }
        if(r.printable!=null)
        {
            this.printable = new String(r.printable);
        }
        else
        {
            this.printable = null;
        }
        this.cliques = new ArrayList<>();
        if(r.cliques!= null)
        {
            for(Clique c: r.cliques)
            {
                this.cliques.add(new Clique(c));
            }
        }
        
    }
    
    /**
     *  Description of the Method
     *
     *@param  file                       Description of the Parameter
     *@exception  IOException            Description of the Exception
     *@exception  FileNotFoundException  Description of the Exception
     */
    public Rcc(String file, String algName, String cliquesFileName) throws IOException 
    {
     int count=0, rcc_class = -1, i=0;
     String line="", command= Parameters.pathTomita+"qc --input-file="+file+" --algorithm="+algName; //"/home/gdelrio/RCCJava/quick-cliques/bin/qc_onebyone --input-file="+file+" --algorithm="+algName;
     String[] data=null;

     BufferedWriter outfileCliques = null;
     boolean flagPrintCliques = false;
     if(cliquesFileName!= null)
     {
        flagPrintCliques = true;
        outfileCliques = new BufferedWriter(new FileWriter(cliquesFileName));
        cliques = new ArrayList<>();
     }

     try
    {
        this.rcc = new ArrayList<>();
        for(i=0; i<26; i++)
         {
             this.rcc.add(0);
         }
        ArrayList<Integer> qc = null;
        Integer v = null;
        Process p = Runtime.getRuntime().exec(command);
        BufferedReader reader1 = new BufferedReader(new InputStreamReader(p.getInputStream()));

        while ((line = reader1.readLine())!= null)
        {
            if(!line.startsWith("NOTE:") && !line.startsWith("Reading") && !line.startsWith("tomita:"))
            {

                data=line.split(" ");
                qc = new ArrayList<>();
                for(i=0; i<data.length; i++)
                {
                    v = new Integer(data[i]);
                    qc.add(v);
                }
                Collections.sort(qc);

                rcc_class = getRCCClass(qc);
              //System.err.println(qc.toString()+" "+rcc_class);

                if(rcc_class>-1)
                {
                  //count = Integer.parseInt((String)this.rcc.get(rcc_class)) + 1;
                  count = this.rcc.get(rcc_class) + 1;
                  //this.rcc.set(rcc_class, String.valueOf(count));
                  this.rcc.set(rcc_class, count);
                    if(flagPrintCliques)
                    {
                        outfileCliques.write(line+"\n");
                        cliques.add( new Clique(qc,rcc_class) );
                    }
                }
            }
        }
        reader1.close();

        //System.out.println(rcc.toString());
        this.printable = rcc.toString();
        if(flagPrintCliques)
        {
            outfileCliques.close();
        }
    }
     catch(Exception e)
      {
       System.err.println("Error GetRCCFromPDB.printRCC:"+e.getMessage());
//           return "Error GetRCCFromPDB.printRCC:"+e.getMessage();
      }
    }
    // end of printRCC method

    /**
     * Description of the Method
     *
     *@param  rcc                       Description of the Parameter
     *@exception  IOException            Description of the Exception
     *@exception  FileNotFoundException  Description of the Exception
     */
/*
    public static int getRCCClass(ArrayList<Integer> rcc)
    {
     int result = -1;
     try
      
     catch(Exception e)
      {
       System.err.println("Error GetRCCFromPDB.getRCCClass:"+e.getMessage());
      }

     return result;
    }
    // end of getRCCClass method
*/
  public static int getRCCClass(ArrayList<Integer> rcc)
  {
    int result = -1;
    try
    {
      int val = 0;int diferencia = 0;int i = 0;int contiguo = 0;
      LinkedList contiguos = new LinkedList();
      for (i = 0; i < rcc.size(); i++)
      {
        if (i == 0)
        {
          val = ((Integer)rcc.get(i)).intValue();
          contiguo = 1;
        }
        else
        {
          diferencia = ((Integer)rcc.get(i)).intValue() - val;
          if (diferencia == 1) { contiguo += 1;
          }
          else {
            if (contiguo > 1) contiguos.add(new Integer(contiguo)); else
              contiguos.add(new Integer(1));
            contiguo = 1;
          }
          val = ((Integer)rcc.get(i)).intValue();
        }
      }
      if (contiguo == 1) contiguos.add(new Integer(1)); else
        contiguos.add(new Integer(contiguo));
      Collections.sort(contiguos);
      


      if (rcc.size() == 3)
      {
        if (((Integer)contiguos.get(0)).intValue() == 3) { result = 2;
        } else if (((Integer)contiguos.get(1)).intValue() == 2) result = 1; else {
          result = 0;
        }
      } else if (rcc.size() == 4)
      {
        if (((Integer)contiguos.get(0)).intValue() == 4) { result = 7;
        } else if (((Integer)contiguos.get(1)).intValue() == 3) { result = 6;
        } else if ((((Integer)contiguos.get(0)).intValue() == 2) && (((Integer)contiguos.get(1)).intValue() == 2)) { result = 5;
        } else if ((((Integer)contiguos.get(0)).intValue() == 1) && (((Integer)contiguos.get(2)).intValue() == 2)) result = 4; else {
          result = 3;
        }
      } else if (rcc.size() == 5)
      {
        if (((Integer)contiguos.get(0)).intValue() == 5) { result = 14;
        } else if ((((Integer)contiguos.get(0)).intValue() == 2) && (((Integer)contiguos.get(1)).intValue() == 3)) { result = 13;
        } else if ((((Integer)contiguos.get(1)).intValue() == 2) && (((Integer)contiguos.get(2)).intValue() == 2)) { result = 12;
        } else if (((Integer)contiguos.get(1)).intValue() == 4) { result = 11;
        } else if ((((Integer)contiguos.get(1)).intValue() == 1) && (((Integer)contiguos.get(2)).intValue() == 3)) { result = 10;
        } else if ((((Integer)contiguos.get(2)).intValue() == 1) && (((Integer)contiguos.get(3)).intValue() == 2)) result = 9; else {
          result = 8;
        }
      } else if (rcc.size() == 6)
      {
        if (((Integer)contiguos.get(0)).intValue() == 6) { result = 25;
        } else if ((((Integer)contiguos.get(0)).intValue() == 3) && (((Integer)contiguos.get(1)).intValue() == 3)) { result = 24;
        } else if ((((Integer)contiguos.get(0)).intValue() == 2) && (((Integer)contiguos.get(1)).intValue() == 4)) { result = 23;
        } else if ((((Integer)contiguos.get(0)).intValue() == 2) && (((Integer)contiguos.get(2)).intValue() == 2)) { result = 22;
        } else if (((Integer)contiguos.get(1)).intValue() == 5) { result = 21;
        } else if ((((Integer)contiguos.get(1)).intValue() == 2) && (((Integer)contiguos.get(2)).intValue() == 3)) { result = 20;
        } else if ((((Integer)contiguos.get(0)).intValue() == 1) && (((Integer)contiguos.get(2)).intValue() == 4)) { result = 19;
        } else if ((((Integer)contiguos.get(2)).intValue() == 2) && (((Integer)contiguos.get(3)).intValue() == 2)) { result = 18;
        } else if ((((Integer)contiguos.get(2)).intValue() == 1) && (((Integer)contiguos.get(3)).intValue() == 3)) { result = 17;
        } else if ((((Integer)contiguos.get(3)).intValue() == 1) && (((Integer)contiguos.get(4)).intValue() == 2)) result = 16; else {
          result = 15;
        }
      }
    }
    catch (Exception e) {
      System.err.println("Error GetRCCFromPDB.getRCCClass:" + e.getMessage());
    }
    
    return result;
  }
      

    public void printTemplates(String outputFileName, int templateLength, String proteinSequence) throws IOException
    {
        BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
        //out.write(this.printable+"\n");
        out.write("sequence,rcc\n");
        
        if( proteinSequence.length() < templateLength )
        {
            out.write(proteinSequence+","+this.printable.replaceAll("\\[", "").replaceAll("\\]", "")+"\n");
            out.close();
            return;
        }
        
        
        int skipSize = templateLength /2;
        int numTemplates = (int) Math.floor((float)proteinSequence.length() / skipSize );
        template[] templates = new template[numTemplates];
        
        int start = 0;
        for(int i=0;i<numTemplates;i++)
        {
            templates[i] = new template();
            templates[i].start = start;
            templates[i].end   = Math.min( start+templateLength, proteinSequence.length() );
            start += skipSize;
        }
        int cNumer =0;
        for(Clique c : cliques)
        {
            if(c.nodes.length < 3)
            {
                continue;
            }
            
            cNumer++;
            int firstTemplate = (int) Math.floor( (float) c.start/skipSize );
            int lastTemplate  = (int) Math.floor( (float) c.end  /skipSize );
            
            //skip cliques spanning more than the templateLength
            if(lastTemplate - firstTemplate > 1)
            {
                continue;
            }
            
            lastTemplate = Math.min(lastTemplate, numTemplates-1);
            
            if(numTemplates == 1)
            {
                lastTemplate = firstTemplate;
            }
            
            for(int t= firstTemplate; t<=lastTemplate; t++)
            {
                try
                {
                    if(templates[t].start <= c.start && templates[t].end >= c.end)
                    {

                        templates[t].rcc[c.rccEntry] = templates[t].rcc[c.rccEntry] + 1;
                    }
                    
                }
                catch(ArrayIndexOutOfBoundsException e)
                    {
                        System.out.print("template " +t+ "  rccEntry " + c.rccEntry + " clique " + cNumer);
                        for(int n : c.nodes)
                        {
                            System.out.print( "  " + n );
                        }
                        System.out.println();
                        throw e;
                    }
            }
        }
        
        start = 0;
        for(int i=0; i<numTemplates;i++)
        {
            int end = Math.min( start+templateLength, proteinSequence.length() );
            templates[i].sequence = proteinSequence.substring( start, end  );
            start+= skipSize;
            out.write(templates[i].sequence +","+ templates[i].getRcc()+"\n");
        }
        out.close();
    }

    public Stats getStats(int numQuantiles)
    {
        ArrayList<Integer> data = new ArrayList<>();
        for(Clique c: cliques)
        {
            data.add(c.end - c.start);
        }
        if(data.size()==0)
        {
            return null;
        }
        
        return new Stats(data, numQuantiles);
    }
    
    public double distanceDif(Rcc r)
    {
        double d = 0;
        for(int i=0; i< 26; i++)
        {
            d += Math.abs(this.rcc.get(i) - r.rcc.get(i));
        }
        return d;
    }
    
    public double distanceDifWeighted(Rcc r)
    {
        double d = 0;
        for(int i=0; i<=2; i++)
        {
            d += Math.abs(this.rcc.get(i) - r.rcc.get(i))*3;
        }
        for(int i=3; i<=7; i++)
        {
            d += Math.abs(this.rcc.get(i) - r.rcc.get(i))*4;
        }
        for(int i=8; i<=14; i++)
        {
            d += Math.abs(this.rcc.get(i) - r.rcc.get(i))*5;
        }
        for(int i=15; i<=25; i++)
        {
            d += Math.abs(this.rcc.get(i) - r.rcc.get(i))*6;
        }
        return d;
    }
}
