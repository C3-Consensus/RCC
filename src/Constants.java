
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.Vector;
import static jdk.nashorn.internal.objects.NativeArray.map;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Constants {
    
    public static Random rand;
    
    public static ArrayList< ArrayList<Double> > anglePriorsPhi;
    public static ArrayList< ArrayList<Double> > anglePriorsPsi;
    public static ArrayList< ArrayList<Double> > anglePriorsOmega;
    public static Map<String, Integer> priorId;
    public static ArrayList<Character> aminoacidList;
    public static ArrayList<String> aminoacidPairList;
    
    public static int numQuantiles = 1000;
    public static double numSearches;
    
    public static double infNeg = -1000;
    public static double infPos = 1000;
    
    public static double distNCa = 1.47;
    public static double distCN  = 1.32;
    public static double distCaC = 1.53;
    public static double distCO = 1.24;
    
    public static double angleDegaCNCa = 180 - 123;//57
    public static double angleDegCaCN = 360 - 180 + 114;//294 ~ 66
    public static double angleDegNCaC = 180 - 110;// 70
    public static double angleDegCaCO = 180 - 121;
    
    public static double angleCNCa = angleDegaCNCa * Math.PI / 180;
    public static double angleCaCN = angleDegCaCN * Math.PI / 180;
    public static double angleNCaC = angleDegNCaC * Math.PI / 180;
    public static double angleCaC0 = angleDegCaCO * Math.PI / 180;
    
    public static double sinAngleCNCa  = Math.sin( angleCNCa * Math.PI / 180);
    public static double sinAngleCaCN = Math.sin( angleCaCN * Math.PI / 180);
    public static double sinAngleNCaC = Math.sin( angleNCaC * Math.PI / 180);
    
    public static double cosAngleCNCa = Math.cos( angleCNCa * Math.PI / 180); 
    public static double cosAngleCaCN = Math.cos(angleCaCN * Math.PI / 180);
    public static double cosAngleNCaC = Math.cos(angleNCaC * Math.PI / 180);           

    public static double angleDegaCNCaFlip = 360 - angleDegaCNCa;
    public static double angleDegCaCNFlip = angleDegCaCN;
    public static double angleDegNCaCFlip = 360 - angleDegNCaC;
    public static double angleDegCaCOFlip = 360 - angleDegCaCO;
    
    public static double angleCNCaFlip = angleDegaCNCaFlip * Math.PI / 180;
    public static double angleCaCNFlip = angleDegCaCNFlip * Math.PI / 180;
    public static double angleNCaCFlip = angleDegNCaCFlip * Math.PI / 180;
    public static double angleCaCOFlip = angleDegCaCOFlip * Math.PI / 180;
    
    public static double sinAngleCNCaFlip  = Math.sin( angleCNCa * Math.PI / 180);
    public static double sinAngleCaCNFlip = Math.sin( angleCaCN * Math.PI / 180);
    public static double sinAngleNCaCFlip = Math.sin( angleNCaC * Math.PI / 180);
    
    public static double cosAngleCNCaFlip = Math.cos( angleCNCa * Math.PI / 180); 
    public static double cosAngleCaCNFlip = Math.cos(angleCaCN * Math.PI / 180);
    public static double cosAngleNCaCFlip = Math.cos(angleNCaC * Math.PI / 180);           
    
    public static final double[] distCaO = new double[26];
    public static final double[] angleDegOCa0Ca1 = new double[26];
    public static final double[] angleCa0OCa1 = new double[26];
    public static final double[] angleDegOCa0Ca1Flip = new double[26];
    public static final double[] angleCa0OCa1Flip = new double[26];
    
    
    static Vector p_charged_aa = new Vector();
    static Vector n_charged_aa = new Vector();
    static Vector hydrophobs_aa = new Vector();
    
    static Map<String,Character> aminoacidSymbol = new HashMap<>();
    
    

    //initializa arrays
    static
    {
        numSearches = Math.log10(numQuantiles)/Math.log10(2);
        
        rand = new Random(1234);
        aminoacidSymbol.put("ALA",'A');
        aminoacidSymbol.put("CYS",'C');
        aminoacidSymbol.put("ASP",'D');
        aminoacidSymbol.put("GLU",'E');
        aminoacidSymbol.put("PHE",'F');
        aminoacidSymbol.put("GLY",'G');
        aminoacidSymbol.put("HIS",'H');
        aminoacidSymbol.put("ILE",'I');
        aminoacidSymbol.put("LYS",'K');
        aminoacidSymbol.put("LEU",'L');
        aminoacidSymbol.put("MET",'M');
        aminoacidSymbol.put("ASN",'N');
        aminoacidSymbol.put("PRO",'P');
        aminoacidSymbol.put("GLN",'Q');
        aminoacidSymbol.put("ARG",'R');
        aminoacidSymbol.put("SER",'S');
        aminoacidSymbol.put("THR",'T');
        aminoacidSymbol.put("VAL",'V');
        aminoacidSymbol.put("TRP",'W');
        aminoacidSymbol.put("TYR",'Y');

        
        p_charged_aa.add("ASP");
        p_charged_aa.add("GLU");
        p_charged_aa.add("ASN");
        p_charged_aa.add("GLN");
        p_charged_aa.add("TYR");
        n_charged_aa.add("ARG");
        n_charged_aa.add("LYS");
        n_charged_aa.add("ASN");
        n_charged_aa.add("GLN");
        n_charged_aa.add("PHE");
        n_charged_aa.add("HIS");
        hydrophobs_aa.add("ALA");
        hydrophobs_aa.add("LEU");
        hydrophobs_aa.add("ILE");
        hydrophobs_aa.add("VAL");
        hydrophobs_aa.add("TRP");
        hydrophobs_aa.add("SER");
        hydrophobs_aa.add("THR");
        hydrophobs_aa.add("CYS");
        hydrophobs_aa.add("MET");

        distCaO[0] = 2.401;  //ala	a
        distCaO[2] = 2.397;  //cys	c
        distCaO[3] = 2.396;  //asp	d
        distCaO[4] = 2.396;  //glu	e
        distCaO[5] = 2.399;  //phe	f
        distCaO[6] = 2.401;  //gly	g
        distCaO[7] = 2.398;  //his	h
        distCaO[8] = 2.4;  //ile	i
        distCaO[10] = 2.4;  //lys	k
        distCaO[11] = 2.395;  //leu	l
        distCaO[12] = 2.399;  //met	m
        distCaO[13] = 2.395;  //asn	n
        distCaO[15] = 2.405;  //pro cis	p
        distCaO[16] = 2.398;  //gln	q
        distCaO[17] = 2.396;  //arg	r
        distCaO[18] = 2.396;  //ser	s
        distCaO[19] = 2.397;  //thr	t
        distCaO[21] = 2.401;  //val	v
        distCaO[22] = 2.403;  //trp	w
        distCaO[23] = 2.397;  //tyr	y
        //distCaO[] = 2.397;  //cyx	
        //distCaO[] = 2.406;  //pro trans	

        angleDegOCa0Ca1[0] = 46.97;  //ala	a
        angleDegOCa0Ca1[2] = 47.07;  //cys	c
        angleDegOCa0Ca1[3] = 47.12;  //asp	d
        angleDegOCa0Ca1[4] = 47.22;  //glu	e
        angleDegOCa0Ca1[5] = 47.89;  //phe	f
        angleDegOCa0Ca1[6] = 47.38;  //gly	g
        angleDegOCa0Ca1[7] = 47.09;  //his	h
        angleDegOCa0Ca1[8] = 42.78;  //ile	i
        angleDegOCa0Ca1[10] = 46.96;  //lys	k
        angleDegOCa0Ca1[11] = 47.22;  //leu	l
        angleDegOCa0Ca1[12] = 47.12;  //met	m
        angleDegOCa0Ca1[13] = 47.2;  //asn	n
        angleDegOCa0Ca1[15] = 46.52;  //pro cis	p
        angleDegOCa0Ca1[16] = 47.27;  //gln	q
        angleDegOCa0Ca1[17] = 47.12;  //arg	r
        angleDegOCa0Ca1[18] = 47.15;  //ser	s
        angleDegOCa0Ca1[19] = 47.11;  //thr	t
        angleDegOCa0Ca1[21] = 47.21;  //val	v
        angleDegOCa0Ca1[22] = 46.91;  //trp	w
        angleDegOCa0Ca1[23] = 47.64;  //tyr	y
        //angleDegOCa0Ca1[] = 47.31;  //cyx	
        //angleDegOCa0Ca1[] = 47.45;  //pro trans	
        
        
    //public static double[] angleCa0OCa1 = new double[26];
    //public static double[] angleDegOCa0Ca1Flip = new double[26];
    //public static double[] angleCa0OCa1Flip = new double[26];
    
        for(int i=0; i<26; i++)
        {
            angleDegOCa0Ca1[i] = angleDegOCa0Ca1[i];
            angleCa0OCa1[i] = angleDegOCa0Ca1[i] * Math.PI / 180;
            angleDegOCa0Ca1Flip[i] = 360 - angleDegOCa0Ca1[i];
            angleCa0OCa1Flip[i] = angleDegOCa0Ca1Flip[i] * Math.PI / 180;
        }
    }

    public static boolean complement(String aa1, String aa2, String pairs) {
            boolean result = false;

            if (pairs.equalsIgnoreCase("all_complementary")) {
                    if ((hydrophobs_aa.contains(aa1) && hydrophobs_aa.contains(aa2)) || (p_charged_aa.contains(aa1) && n_charged_aa.contains(aa2)) || (p_charged_aa.contains(aa2) && n_charged_aa.contains(aa1))) {
                            result = true;
                    }
            } else if (pairs.equalsIgnoreCase("h_complementary")) {
//			System.err.println(pairs + ":" + aa1 + " " + aa2);

                    if (hydrophobs_aa.contains(aa1) && hydrophobs_aa.contains(aa2)) {
//				System.err.println("got it");
                            result = true;
                    }
            } else if (pairs.equalsIgnoreCase("ch_complementary")) {
                    if ((p_charged_aa.contains(aa1) && n_charged_aa.contains(aa2)) || (p_charged_aa.contains(aa2) && n_charged_aa.contains(aa1))) {
                            result = true;
                    }
            } else {
                    return false;
            }
            return result;
    }
    
    public static int aminoacidPairToId(String aa)
    {
        aa = aa.toLowerCase();
        int i = -1, j=-1;
        
        switch(aa.charAt(0))
        {
            case 'a':
                i = 0;
                break;
            case 'c':
                i = 1;
                break;
            case 'd':
                i = 2;
                break;
            case 'e':
                i = 3;
                break;
            case 'f':
                i = 4;
                break;
            case 'g':
                i = 5;
                break;
            case 'h':
                i = 6;
                break;
            case 'i':
                i = 7;
                break;
            case 'k':
                i = 8;
                break;
            case 'l':
                i = 9;
                break;
            case 'm':
                i = 10;
                break;
            case 'n':
                i = 11;
                break;
            case 'p':
                i = 12;
                break;
            case 'q':
                i = 13;
                break;
            case 'r':
                i = 14;
                break;
            case 's':
                i = 15;
                break;
            case 't':
                i = 16;
                break;
            case 'v':
                i = 17;
                break;
            case 'w':
                i = 18;
                break;
            case 'y':
                i = 19;
                break;
        }
        switch(aa.charAt(1))
        {
            case 'a':
                j = 0;
                break;
            case 'c':
                j = 1;
                break;
            case 'd':
                j = 2;
                break;
            case 'e':
                j = 3;
                break;
            case 'f':
                j = 4;
                break;
            case 'g':
                j = 5;
                break;
            case 'h':
                j = 6;
                break;
            case 'i':
                j = 7;
                break;
            case 'k':
                j = 8;
                break;
            case 'l':
                j = 9;
                break;
            case 'm':
                j = 10;
                break;
            case 'n':
                j = 11;
                break;
            case 'p':
                j = 12;
                break;
            case 'q':
                j = 13;
                break;
            case 'r':
                j = 14;
                break;
            case 's':
                j = 15;
                break;
            case 't':
                j = 16;
                break;
            case 'v':
                j = 17;
                break;
            case 'w':
                j = 18;
                break;
            case 'y':
                j = 19;
                break;
        }
        return j*20+i;
        
    }
    
    public static void initializePriors(String path) throws FileNotFoundException, IOException
    {
        anglePriorsOmega = new ArrayList<>();
        anglePriorsPhi = new ArrayList<>();
        anglePriorsPsi = new ArrayList<>();
        priorId = new HashMap<>();

        aminoacidPairList = new ArrayList<>();
        aminoacidList = new ArrayList<>();
        aminoacidList.add('a');
        aminoacidList.add('c');
        aminoacidList.add('d');
        aminoacidList.add('e');
        aminoacidList.add('f');
        aminoacidList.add('g');
        aminoacidList.add('h');
        aminoacidList.add('i');
        aminoacidList.add('k');
        aminoacidList.add('l');
        aminoacidList.add('m');
        aminoacidList.add('n');
        aminoacidList.add('p');
        aminoacidList.add('q');
        aminoacidList.add('r');
        aminoacidList.add('s');
        aminoacidList.add('t');
        aminoacidList.add('v');
        aminoacidList.add('w');
        aminoacidList.add('y');

        
        for(int i=0; i<aminoacidList.size(); i++)
        {
            for(int j=0; j<aminoacidList.size(); j++)
            {
                aminoacidPairList.add(""+aminoacidList.get(i)+aminoacidList.get(j));
                priorId.put(""+ aminoacidList.get(j)+aminoacidList.get(i), j*20+i );
            }
        }
        String line;
        ArrayList<Double> prior;
        BufferedReader in;
        for(String aa : aminoacidPairList)
        {
            //PHI
            prior = new ArrayList<>();
            in = new BufferedReader(new FileReader(path+"phi/quantiles/ "+aa+" .csv"));
            in.readLine();
            while( (line = in.readLine()) != null)
            {
                prior.add( Double.parseDouble( line.split(",")[1] )  );
            }
            in.close();
            anglePriorsPhi.add(prior);

            //PSI
            prior = new ArrayList<>();
            in = new BufferedReader(new FileReader(path+"psi/quantiles/ "+aa+" .csv"));
            in.readLine();
            while( (line = in.readLine()) != null)
            {
                prior.add( Double.parseDouble( line.split(",")[1] )  );
            }
            in.close();
            anglePriorsPsi.add(prior);

            //OMEGA
            prior = new ArrayList<>();
            in = new BufferedReader(new FileReader(path+"omega/quantiles/ "+aa+" .csv"));
            in.readLine();
            while( (line = in.readLine()) != null)
            {
                prior.add( Double.parseDouble( line.split(",")[1] )  );
            }
            in.close();
            anglePriorsOmega.add(prior);            
        }
        
        
    }
    
}
