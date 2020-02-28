
import java.util.ArrayList;
import java.util.Collections;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class Parameters {
    public static double min   = 0;
    public static double max   = 5;
    public static String distance_criterium = "once";
    public static String pairs = "all";
    public static Boolean flagLateralChain = false;
    public static String pathTomita = "";
    public static String graphFileName = "";
    public static String cliquesFileName = "";
    public static String pdbFileName = "";
    public static boolean loadCorrect = false;
    public static String chain = "";
    public static String fileList = null;
    public static String pathPdb = null;
    
    public static void LoadParameters(String[] args)
    {
        //command type: uniprot, cath domain, annealing
        
        ArrayList<String> options = new ArrayList<>();
        Collections.addAll(options, args);
        int optionIndex;
        
        Parameters.loadCorrect = true;
        //display options
        if(args.length==0)
        {
            System.out.println("-pdb         \t<file to be processed>");
            System.out.println("-chain     \t<chain in the pdb file to be processed>");
            System.out.println("-pathTomita\t<path to tomita executable file>");
            System.out.println("-lateralChain  \t<yes/no (no)>");
            System.out.println("-minDist   \t<minimum contact distance (0)>");
            System.out.println("-maxDist   \t<maximum contact distance (7)>");
            System.out.println("-fileList  \t<csv file with filename, chain>");
            System.out.println("-pathPdbFolder   \t<path to pdb database folder>");
            
            Parameters.loadCorrect = false;
            return;
        }

        optionIndex = options.indexOf("-fileList")+1;
        if(optionIndex ==0)
        {
            System.out.println("no fileList");      
        }
        else
        {
            Parameters.fileList = options.get(optionIndex);
            optionIndex = options.indexOf("-pathPdbFolder")+1;
            if(optionIndex == 0)
            {
                System.out.println("-pathPdbFolder missing");
                Parameters.loadCorrect = false;
                return;
            }
            Parameters.pathPdb = options.get(optionIndex);
        }

        
        //System.out.println("-pathTomita path to tomita executable file");
            optionIndex = options.indexOf("-pdb") + 1;
            if(optionIndex == 0 && Parameters.fileList == null)
            {
                System.out.println("-pdb missing");
                Parameters.loadCorrect = false;
                return;
            }
            else
            {
                Parameters.pdbFileName = options.get(optionIndex);
            }
            
            
            //System.out.println("-pathTomita path to tomita executable file");
            optionIndex = options.indexOf("-chain") + 1;
            if(optionIndex == 0 && Parameters.fileList == null)
            {
                System.out.println("-chain missing");
                Parameters.loadCorrect = false;
                return;
            }
            else
            {
                Parameters.chain = options.get(optionIndex);
            }
            
            //System.out.println("-pathTomita path to tomita executable file");
            optionIndex = options.indexOf("-pathTomita") + 1;
            if(optionIndex == 0)
            {
                System.out.println("-pathTomita missing");
                Parameters.loadCorrect = false;
                return;
            }
            Parameters.pathTomita = options.get(optionIndex);
            

            //System.out.println("-backbone yes/no (yes)");
            optionIndex = options.indexOf("-backbone") + 1;
            if(optionIndex == 0)
            {
                Parameters.flagLateralChain = false;
            }
            else
            {
                Parameters.flagLateralChain = options.get(optionIndex).equalsIgnoreCase("yes");
            }
            
            //System.out.println("-minDist minimum contact distance (0)");
            optionIndex = options.indexOf("-minDist") + 1;
            if(optionIndex == 0)
            {
                Parameters.min = 0;
            }
            else
            {
                Parameters.min = Integer.parseInt(options.get(optionIndex));
            }
            
            //System.out.println("-maxDist maximum contact distance (5)");
            optionIndex = options.indexOf("-maxDist") + 1;
            if(optionIndex == 0)
            {
                Parameters.max = 7;
            }
            else
            {
                Parameters.max = Integer.parseInt(options.get(optionIndex));
            }
            
            Parameters.graphFileName = Parameters.pdbFileName + "graph";
            
            System.out.println("Path to tomita binary: " + Parameters.pathTomita);
            System.out.println("Using lateral chain:" + Parameters.flagLateralChain);
            System.out.println("Contact minimum distance: " + Parameters.min);
            System.out.println("Contact maximum distance: " + Parameters.max);
            System.out.println("Input pdb file: " + Parameters.pdbFileName);
            System.out.println("File list: " + Parameters.fileList);
            System.out.println("Path to pdb folder " + Parameters.pathPdb);
    }
}
