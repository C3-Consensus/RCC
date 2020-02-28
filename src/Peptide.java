import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Peptide 
{
    String chain;
    Reticle reticle;
    ArrayList<Residue> residues;
    ArrayList<Integer> residuesInitialized;
    String sequence;
    
    //not used?
    ArrayList<String> atoms;
    
    ArrayList<Angles> angles;
    String comments;
    int missingResidues;
    Point3d upperBound;
    Point3d lowerBound;
    
    //generate empty peptide
    public Peptide()
    {
        reticle = null;
        residues = new ArrayList<>();
        residuesInitialized = new ArrayList<>();
        atoms = new ArrayList<>();
        angles = new ArrayList<>();
        comments = "";
        missingResidues = 0;
        sequence = null;
    }
    
    public Peptide(Peptide p)
    {
        if(p.chain != null)
        {
            chain = new String(p.chain);
        }
        else
        {
            chain = null;
        }
        if(p.sequence != null)
        {
            sequence = new String(p.sequence);            
        }
        else
        {
            sequence = null;
        }

        if(p.reticle != null)
        {
            reticle = new Reticle(p.reticle);
        }
        else
        {
            reticle = null;
        }
        if(p.residues != null)
        {
            residues = new ArrayList<>();
            for(Residue r : p.residues)
            {
                residues.add(new Residue(r));
            }
        }
        if(p.residuesInitialized != null)
        {
            residuesInitialized = new ArrayList<>();
            for(Integer i : p.residuesInitialized)
            {
                residuesInitialized.add(new Integer(i));
            }
        }
        else
        {
            residuesInitialized = null;
        }
        
        if(p.atoms != null)
        {
            atoms = new ArrayList<>();
            for(String at: p.atoms)
            {
                atoms.add(new String(at));
            }
        }
        if(p.angles != null)
        {
            angles = new ArrayList<>();
            for(Angles t : p.angles)
            {
                angles.add(new Angles(t));
            }
        }
        if(p.comments != null)
        {
            comments = new String(p.comments);
        }
        else
        {
            comments = null;
        }
        missingResidues = p.missingResidues;
        if(p.upperBound != null)
        {
            upperBound = new Point3d(p.upperBound);
        }
        else
        {
            upperBound = null;
        }
        if(p.lowerBound != null)
        {
            lowerBound = new Point3d(p.lowerBound);
        }
        else
        {
            lowerBound = null;
        }
    }
    
    public Peptide(ArrayList<Angles> angles)
    {
        reticle = null;
        residues = new ArrayList<>();
        residuesInitialized = null;
        atoms = null;
        this.angles = new ArrayList<>();
        for(Angles a : angles)
        {
            this.angles.add(new Angles(a));
        }
        chain = null;
        comments = null;
        
    }
    
    
    //generate peptide from set of dihedral angles
    public Peptide(ArrayList<Angles> angles, ArrayList<Character> aminoacids)
    {
        reticle = null;
        residues = new ArrayList<>();
        residuesInitialized = new ArrayList<>();
        atoms = new ArrayList<>();
        angles = new ArrayList<>();
        comments = "";
        
        this.angles = angles;

        Residue r1 = new Residue(1, aminoacids.get(0));

        r1.atomsBackbone.N = new Atom(1,"N");
        r1.atomsBackbone.Ca = new Atom(2,"CA");
        r1.atomsBackbone.C = new Atom(3,"C");
        
        
        r1.atomsBackbone.N.coords.setCoords(0, 0, 0);
        r1.atomsBackbone.N.name = "N";
        
        r1.atomsBackbone.Ca.coords.setCoords(Constants.distNCa, 0, 0);
        r1.atomsBackbone.Ca.name = "Ca";
        
//        r1.atomsBackbone.C.coords.setCoords( Constants.cosAngleNCaC, Constants.sinAngleNCaC, 0);
        
        r1.atomsBackbone.C.coords.setCoords( 1, 0, 0);
        r1.atomsBackbone.C.coords.rotateWithRodriguez(Constants.angleNCaC, new Point3d(0,0,-1));
        
        r1.atomsBackbone.C.coords.multiplyByScalar(Constants.distCaC);
        r1.atomsBackbone.C.coords.sumVector(r1.atomsBackbone.Ca.coords);
        r1.atomsBackbone.C.name = "C";
         
        r1.dihedralAngles = angles.get(0);
        
        r1.vectorNCa = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.Ca.coords, r1.atomsBackbone.N.coords);
        r1.vectorNCa.normalize();
        r1.vectorCaC = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.C.coords, r1.atomsBackbone.Ca.coords);
        r1.vectorCaC.normalize();

        this.residues.add(r1);

        boolean flipFlag = false;
        Residue prev = r1;

        for(int r = 1; r< angles.size(); r++)
        {
            flipFlag = !flipFlag;
            Residue res = new Residue(r, aminoacids.get(r));
            res.setResidue(prev, angles.get(r), flipFlag);
            prev = res;
            this.residues.add(res);
        }
    }
    
    //get peptide from pdb file with specified chain and calculate reticle containing it. File may be .pdb or .pdb.gz
    public Peptide(String fileName, String chain, boolean flagBackbone) throws IOException, FileNotFoundException, IllegalArgumentException, IllegalAccessException 
    {
        this.reticle = null;
        this.residues = new ArrayList<>();
        this.residuesInitialized = new ArrayList<>();
        this.atoms = new ArrayList<>();
        this.angles = new ArrayList<>();
        this.comments = "";
        this.missingResidues = 0;

        this.chain = chain;
        String line = "";
        String atomName = "";
        String resName = "";
        String chainID = "";
        Integer resSeq = 0;
        Integer serial = 0;
        String pos = "";
        double xS = 0;
        double yS = 0;
        double zS = 0;
        String aa = "";
        String aa_old = "";
        double avgX = 0.0;
        double avgY = 0.0;
        double avgZ = 0.0;
        double minX = 1000;
        double maxX = -1000;
        double minY = 1000;
        double maxY = -1000;
        double minZ = 1000;
        double maxZ = -1000;

        BufferedReader infile;
        
        if(fileName.contains(".gz"))
        {
            //GZIPInputStream gzip = new GZIPInputStream(new FileInputStream("F:/gawiki-20090614-stub-meta-history.xml.gz"));
            infile = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(fileName)) ));
            //infile = new BufferedReader(new Gz FileReader(fileName));            
        }
        else
        {
            infile = new BufferedReader(new FileReader(fileName));
        }
        
        
        int i = 0, count=-1;
        while ((line = infile.readLine()) != null) 
        {
            if (!line.startsWith("ATOM")) 
            {
                continue;
            }
            if( fileName.contains(".cif"))
            {
                while(line.contains("  "))
                {
                    line = line.replaceAll("  ", " ");
                }
                String[] fields = line.split(" ");

                /*
            0     _atom_site.group_PDB 
            1     _atom_site.id 
            2     _atom_site.type_symbol 
            3     _atom_site.label_atom_id 
            4     _atom_site.label_alt_id 
            5     _atom_site.label_comp_id 
            6     _atom_site.label_asym_id 
            7     _atom_site.label_entity_id 
            8     _atom_site.label_seq_id 
            9     _atom_site.pdbx_PDB_ins_code 
            10    _atom_site.Cartn_x 
            11    _atom_site.Cartn_y 
            12    _atom_site.Cartn_z 
            13    _atom_site.occupancy 
            14    _atom_site.B_iso_or_equiv 
            15    _atom_site.pdbx_formal_charge 
            16    _atom_site.auth_seq_id 
            17    _atom_site.auth_comp_id 
            18    _atom_site.auth_asym_id 
            19    _atom_site.auth_atom_id 
            20    _atom_site.pdbx_PDB_model_num 

            0      1      2  3     4 5   6  7  8   9 10      11      12      13   14    15 16  17  18 19    20
            ATOM   86819  N  N     . THR OB 66 1   ? 240.774 299.008 273.887 1.00 85.41  ? 2   THR AJ N     1 
            ATOM   86820  C  CA    . THR OB 66 1   ? 239.733 299.148 272.880 1.00 82.21  ? 2   THR AJ CA    1 
            ATOM   86821  C  C     . THR OB 66 1   ? 240.312 299.572 271.536 1.00 79.25  ? 2   THR AJ C     1 
            ATOM   86822  O  O     . THR OB 66 1   ? 240.452 300.765 271.268 1.00 79.09  ? 2   THR AJ O     1 
            ATOM   86823  C  CB    . THR OB 66 1   ? 238.670 300.167 273.310 1.00 20.00  ? 2   THR AJ CB    1 
                */
                atomName = fields[3];

                if( !flagBackbone )
                {
                    // skip non-backbone atoms
                    if( !atomName.equals("C") && !atomName.equals("O") && !atomName.equals("N") && !atomName.equals("CA"))
                    {
                        continue;
                    }
                }

                    //aminoacid name: PRO, GLN, VAL, ...
                resName = fields[5];
                    //residue number 
                resSeq = Integer.parseInt(fields[8]);
                    //atom number
                /*
                if(flagBackbone)
                 {
                  if(atomName.equals("C") || atomName.equals("O") || atomName.equals("N") || atomName.equals("CA"))
                     serial = Integer.parseInt(line.substring(6, 11).trim());
                 }
                else serial = Integer.parseInt(line.substring(6, 11).trim());
                */
                serial = Integer.parseInt(fields[1]);
                xS = Double.parseDouble(fields[10]);
                yS = Double.parseDouble(fields[11]);
                zS = Double.parseDouble(fields[12]);

                chainID = fields[18];
                
                if(!(chainID.equalsIgnoreCase(chain) || chain.equalsIgnoreCase("none") || chain.equalsIgnoreCase("all")) ) 
                {
                    continue;
                }


            }
            else
            {
                chainID = line.substring(21, 22);
                if(!(chainID.equalsIgnoreCase(chain) || chain.equalsIgnoreCase("none") || chain.equalsIgnoreCase("all")) ) 
                {
                    continue;
                }

//           1         2         3         4         5         6         7         8
//  1234567890123456789012345678901234567890123456789012345678901234567890123456789012                                    
//  ATOM      1  N   PRO A   1      -3.190   7.728  33.820  1.00 21.66           N
//  ATOM     15  OE1 GLN A   2       5.366   6.856  36.946  1.00 31.80           O  
                            //atom name: N, CA, C, O, CB, CG, CD, ...
                atomName = line.substring(12, 16).trim();
//IGNORE ALTERNATIVE POSITION OF ATOMS!!!
                if(line.charAt(16)!= ' ' && line.charAt(16)!= 'A')
                {
                    continue;
                }
                if( !flagBackbone )
                {
                    // skip non-backbone atoms
                    if( !atomName.equals("C") && !atomName.equals("O") && !atomName.equals("N") && !atomName.equals("CA"))
                    {
                        continue;
                    }
                }

                    //aminoacid name: PRO, GLN, VAL, ...
                resName = line.substring(17, 20);
                    //residue number 
                resSeq = Integer.parseInt(line.substring(22, 26).trim());
                    //atom number
                serial=0;
                /*
                if(flagBackbone)
                 {
                  if(atomName.equals("C") || atomName.equals("O") || atomName.equals("N") || atomName.equals("CA"))
                     serial = Integer.parseInt(line.substring(6, 11).trim());
                 }
                else serial = Integer.parseInt(line.substring(6, 11).trim());
                */
                serial = Integer.parseInt(line.substring(6, 11).trim());
                xS = Double.parseDouble(line.substring(30, 38).trim());
                yS = Double.parseDouble(line.substring(38, 46).trim());
                zS = Double.parseDouble(line.substring(46, 54).trim());
            }
            //get bounds for the reticule
            if( xS < minX )
            {
                minX = xS;
            }
            if( xS > maxX )
            {
                maxX = xS;
            }
            if( yS < minY )
            {
                minY = yS;
            }
            if( yS > maxY )
            {
                maxY = yS;
            }
            if( zS < minZ )
            {
                minZ = zS;
            }
            if( zS > maxZ )
            {
                maxZ = zS;
            }

            if(!serial.equals(""))
            {
                //search can be skipped if the .pdb file is sorted by residue id
                Residue currentResidue;
                int residueId = this.residuesInitialized.indexOf(resSeq);
                if(residueId == -1)
                {
                    this.residuesInitialized.add(resSeq);
                    currentResidue = new Residue(resSeq, resName);
                    //currentResidue.id = resSeq;
                    //currentResidue.aminoacid = resName;
                    //currentResidue.atoms = new ArrayList<>();
                    this.residues.add(currentResidue);
                }
                else
                {
                    currentResidue = this.residues.get(residueId);
                }
                //add atom to residue
                Atom a = new Atom(serial, atomName, xS, yS, zS);
                if(atomName.contentEquals("N") )
                {
                    currentResidue.atomsBackbone.N = a;
                }
                else if(atomName.contentEquals("C") )
                {
                    currentResidue.atomsBackbone.C = a;
                }
                else if(atomName.contentEquals("CA") )
                {
                    currentResidue.atomsBackbone.Ca = a;
                }
                else if(atomName.contentEquals("O"))
                {
                    currentResidue.atomsBackbone.O = a;
                }
                else
                {
                    currentResidue.atomsLateralChain.add(a);
                }
            }
        }

        double maxRadius = 0;

        //calculate center and radius for residues, and asigned reticle member
        for(int r=0; r < this.residues.size(); r++)
        {
            this.residues.get(r).calculateCenterRadiusAndBounds();
            if( maxRadius < this.residues.get(r).radius)
            {
                maxRadius = this.residues.get(r).radius;
            }
        }
        infile.close();

        this.reticle = new Reticle(minX, maxX, minY, maxY, minZ, maxZ);
        this.reticle.maxRadius = maxRadius;
    }
    
        public Graph buildGraph(double min, double max, String distance_criterium, String pairs, boolean flagBackbone, String outname) throws IOException, FileNotFoundException {
		int order = 0;
                
		double xx = 0.0;
		double yy = 0.0;
		double zz = 0.0;
                double distance = 0.0;
                double distanceUpper = 0.0;
                double distanceLower = 0.0;
                
                Graph g = new Graph();
                                
                this.assignReticleToResidues(max);

		if (distance_criterium.equalsIgnoreCase("once"))
                {
                //get residue distance
                    for(int i=0; i< this.residues.size(); i++)
                    {

                    //Members of the same cube in the reticle are automatically close
                    /*    
        //MAY NEED TO REMOVE IF MEMORY IS NOT ENOUGH FOR RETICLE RESOLUTION (DISTANCE FOR RETICLE IS GREATER THAN MAX DISTANCE)
                        for( int j= 0; j< Reticle.cubeMembers.get( residues.get(i).reticleCube ).size(); j++)
                        {
                            //if residues are incompatible, there is no need to calculate distance
                            if( !pairs.equalsIgnoreCase("all") &&  !complement( residues.get(i).aminoacid, residues.get(j).aminoacid, pairs)  )
                            {
                                continue;
                            }
                            
                            //condition to ensure no repeated edges are added
                            if(i < j)
                            {                            
                                edge e = new edge();
                                e.start = i;
                                e.end = j;
                                edges.add(e);
                            }
                        }
                      */  
                    //get distance to neighbours
                        ArrayList<Integer> neighbours = this.reticle.getNeighbours(this.residues.get(i).reticleCube );
                
                        for(int neighbour = 0 ; neighbour < neighbours.size(); neighbour++)
                        {
                            int j = neighbours.get(neighbour);
                                                                                    
                        //skip neighbours that have already been calculated
                            if (j <= i)
                            {
                                continue;
                            }
/*                            
                            if(i == 126 && j==132)
                            {
                                j=j;
                            }
*/
                        //continguous members of the chain are always close
/*                            if( j == i+1)
                            {
                                edge e = new edge();
                                e.start = i;
                                e.end = j;
                                g.edges.add(e);
                                continue;
                            }
  */                          
                            //if residues are incompatible, there is no need to calculate distance
                            if( !pairs.equalsIgnoreCase("all") &&  !Constants.complement( this.residues.get(i).aminoacid.name, this.residues.get(j).aminoacid.name, pairs)  )
                            {
                                continue;
                            }
                            distance = AlgebraFunctions.calculateDistance(this.residues.get(i).coords.x, this.residues.get(i).coords.y, this.residues.get(i).coords.z,
                                                                   this.residues.get(j).coords.x, this.residues.get(j).coords.y, this.residues.get(j).coords.z);

                            //adjust for radius
                            distanceLower = distance - this.residues.get(i).radius - this.residues.get(j).radius;
                            //may be refined to the minimum of adding just one radius
                            distanceUpper = Math.min( distance + this.residues.get(i).radius, distance + this.residues.get(j).radius);

                            //if upper distance is lower than max, all the atoms have to be close
                            if(distanceUpper < max)
                            {
                                Edge e = new Edge();
                                e.start = i;
                                e.end = j;
                                g.edges.add(e);
                                continue;
                            }
                            //if lower distance is higher than max, then no atoms can be close
                            if(distanceLower > max)
                            {
                                continue;
                            }
                            
//                            residue ri = p.residues.get(i);
//                            residue rj = p.residues.get(j);
                            
        //TEST
                            /**/
                    if( distance < Math.max( max + this.residues.get(i).radius, max + this.residues.get(j).radius ) )
                    {
                        double volumeIntersection1 = AlgebraFunctions.volumeIntersection(distance, max, this.residues.get(j).radius);
                        double volumeR2 = 4.0 * Math.PI * Math.pow( this.residues.get(j).radius, 3) / 3.0;
                        
                        double volumeIntersection2 = AlgebraFunctions.volumeIntersection(distance, this.residues.get(i).radius, max);
                        double volumeR1 = 4.0 * Math.PI * Math.pow( this.residues.get(i).radius, 3) / 3.0;


                        double overlapFactor = 0.5;
                        if( volumeIntersection1 > volumeR2 * overlapFactor && volumeIntersection2   > volumeR1 * overlapFactor)
                        {
                            Edge e = new Edge();
                            e.start = i;
                            e.end = j;
                            g.edges.add(e);
                            continue;
                        }
                        
                    }
                            /**/
                            
                            
                            //calculate accurate distance between atoms
                            // far is used as flag to break out of the loop when a pair of atoms that are close are found
                            boolean flagClose = false;
                            boolean flagOnlyBackbone = false;
                            
                            Iterator<Atom> itAtomI = this.residues.get(i).iterator(flagOnlyBackbone);
                            Atom currentAtomI;
                            while( itAtomI.hasNext() && !flagClose )
                            //for(int ia = 0; ia < p.residues.get(i).atoms.size() && !flagClose; ia++)
                            {
                                
                            
                                
                                currentAtomI = itAtomI.next();
                                //distance = Functions.calculateDistance(p.residues.get(i).atoms.get(ia).coords, p.residues.get(j).coords);

                                distance = AlgebraFunctions.calculateDistance( currentAtomI.coords, this.residues.get(j).coords);
                                distanceLower = distance - this.residues.get(j).radius;
                                distanceUpper = distance + this.residues.get(j).radius;
                                
                                //if upper distance is lower than max, all atoms in j have to be close to atom ia of residue i
                                if(distanceUpper < max)
                                {
                                    Edge e = new Edge();
                                    e.start = i;
                                    e.end = j;
                                    g.edges.add(e);
                                    flagClose = true;
                                    continue;
                                }
                                
                                //if lower distance is higher than max, then no atoms in residue j can be close to atom ia of residue i
                                if (distance - this.residues.get(j).radius > max)
                                {
                                    continue;
                                }
                                
                                Iterator<Atom> itAtomJ = this.residues.get(j).iterator(flagBackbone);
                                Atom currentAtomJ;
                                
                                while( itAtomJ.hasNext() && !flagClose)
                                //for(int ja = 0; ja < p.residues.get(j).atoms.size() && !flagClose; ja++)
                                {
                                    currentAtomJ = itAtomJ.next();
                                    distance = AlgebraFunctions.calculateDistance(currentAtomI.coords, currentAtomJ.coords);
                                    //distance = Functions.calculateDistance(p.residues.get(i).atoms.get(ia).coords, p.residues.get(j).atoms.get(ja).coords);
                                    if(distance < max)
                                    {
                                        Edge e = new Edge();
                                        e.start = i;
                                        e.end = j;
                                        g.edges.add(e);
                                        flagClose = true;
                                    }
                                }
                            }
                        }
                    }
		 }

                
                //PENDING
                // if(distance_criterium.equalsIgnoreCase("average"))
                
            g.printGraph(outname, this.residues.size());
            return g;
        }
	// end of doGraph method
         

    //extract domain from current peptide
    public Peptide ExtractDomain(Domain dom)
    {
        Peptide p = new Peptide();
        p.chain = this.chain;
        p.reticle = this.reticle;
        for(int segment = 0; segment < dom.numSegments; segment++)
        {
            for(int i = dom.start[segment]; i <= dom.end[segment]; i++)
            {
                int residueId = this.residuesInitialized.indexOf(i);
            //report residues in domain not found on the pdb file
                if(residueId == -1)
                {
                    p.missingResidues++;
//                    System.err.println("residue not found segment " + segment + " residue " + i );
                    if(!p.comments.equalsIgnoreCase(""))
                    {
                        p.comments+="|";
                    }
                    p.comments += i;
                }
                else
                {
                    p.residuesInitialized.add(i);
                    p.residues.add( this.residues.get(residueId) );
                }
            }
        }
        return p;
    }
    
    //print internal (non-dihedral) angles of the structure
    public void printInternalAngles(String fileName) throws IOException
    {
        Point3d NCa, CaC, CN0, CN1, CO, OC, OCa0, Ca0Ca1;
        Point3d nNCaC, nCaCN, nCNCa, nCaCO, nOCN, nOCa0Ca1;
        double aNCaC,  aCaCN, aCNCa, aCaCO, aOCN, aOCa0Ca1;
        double distCO;

        BufferedWriter outfile = new BufferedWriter(new FileWriter(fileName));
        outfile.write("residue1,residue2,CaCN,CNCa,NCaC,CaCO,OCN,OCa0Ca1,distCO\n");
        
        distCO = AlgebraFunctions.calculateDistance(residues.get(0).atomsBackbone.C.coords, residues.get(0).atomsBackbone.O.coords);
        
        NCa = AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.Ca.coords ,  residues.get(0).atomsBackbone.N.coords );
        CaC = AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.C.coords ,  residues.get(0).atomsBackbone.Ca.coords );
        CO =  AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.O.coords ,  residues.get(0).atomsBackbone.C.coords );
        OC =  AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.C.coords ,  residues.get(0).atomsBackbone.O.coords );
        CN1 = AlgebraFunctions.SubstractPoints( residues.get(1).atomsBackbone.N.coords, residues.get(0).atomsBackbone.C.coords);
        OCa0 = AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.Ca.coords, residues.get(0).atomsBackbone.O.coords);
        Ca0Ca1 = AlgebraFunctions.SubstractPoints( residues.get(1).atomsBackbone.Ca.coords, residues.get(0).atomsBackbone.Ca.coords);
        
        NCa.normalize();
        CaC.normalize();
        CO.normalize();
        OC.normalize();
        CN1.normalize();
        OCa0.normalize();
        Ca0Ca1.normalize();

        nNCaC = AlgebraFunctions.CrossProduct(NCa.getInverse(),CaC);
        nCaCO = AlgebraFunctions.CrossProduct(CaC.getInverse(),CO);
        nOCN = AlgebraFunctions.CrossProduct(OC.getInverse(), CN1);
        nOCa0Ca1 = AlgebraFunctions.CrossProduct(OCa0.getInverse(), Ca0Ca1);
        
        nNCaC.normalize();
        nCaCO.normalize();
        
        double dummy = Math.acos(AlgebraFunctions.InnerProduct(nNCaC, nCaCO));
        
        nOCN.normalize();
        nOCa0Ca1.normalize();
        
        aNCaC = AlgebraFunctions.CalculateAngle(NCa.getInverse(), CaC, nNCaC)*180/Math.PI;
        aCaCO = AlgebraFunctions.CalculateAngle(CaC.getInverse(), CO,  nCaCO)*180/Math.PI;
        aOCN  = AlgebraFunctions.CalculateAngle(OC.getInverse(), CN1, nOCN)*180/Math.PI;
        aOCa0Ca1 = AlgebraFunctions.CalculateAngle(OCa0.getInverse(), Ca0Ca1, nOCa0Ca1)*180/Math.PI;
//        aNCaC = Math.acos( Functions.InnerProduct(NCa, CaC))*180/Math.PI;
        
        
        outfile.write("0,1,NA,NA,"+aNCaC+","+aCaCO+","+aOCN+ "," + aOCa0Ca1 + "," + distCO +"\n");
        
        
        for(int r=1; r<residues.size()-2;r++)
        {
            distCO = AlgebraFunctions.calculateDistance(residues.get(r).atomsBackbone.C.coords, residues.get(r).atomsBackbone.O.coords);
            NCa = AlgebraFunctions.SubstractPoints( residues.get(r).atomsBackbone.Ca.coords, residues.get(r).atomsBackbone.N.coords);
            CaC = AlgebraFunctions.SubstractPoints( residues.get(r).atomsBackbone.C.coords, residues.get(r).atomsBackbone.Ca.coords);
            CO =  AlgebraFunctions.SubstractPoints( residues.get(r).atomsBackbone.O.coords ,  residues.get(r).atomsBackbone.C.coords );
            OC =  AlgebraFunctions.SubstractPoints( residues.get(r).atomsBackbone.C.coords ,  residues.get(r).atomsBackbone.O.coords );
            CN1 = AlgebraFunctions.SubstractPoints( residues.get(r+1).atomsBackbone.N.coords, residues.get(r).atomsBackbone.C.coords);
            OCa0 = AlgebraFunctions.SubstractPoints( residues.get(r).atomsBackbone.Ca.coords, residues.get(r).atomsBackbone.O.coords);
            Ca0Ca1 = AlgebraFunctions.SubstractPoints( residues.get(r+1).atomsBackbone.Ca.coords, residues.get(r).atomsBackbone.Ca.coords);
            CN0 = AlgebraFunctions.SubstractPoints( residues.get(r).atomsBackbone.N.coords, residues.get(r-1).atomsBackbone.C.coords);
        
            nCaCN = AlgebraFunctions.CrossProduct(CaC.getInverse(), CN1);
            nCNCa = AlgebraFunctions.CrossProduct(CN0.getInverse(), NCa);
            nNCaC = AlgebraFunctions.CrossProduct(NCa.getInverse(), CaC);
            nCaCO = AlgebraFunctions.CrossProduct(CaC.getInverse(),CO);
            nOCN = AlgebraFunctions.CrossProduct(OC.getInverse(), CN1);
            nOCa0Ca1 = AlgebraFunctions.CrossProduct(OCa0.getInverse(), Ca0Ca1);
        
            CN0.normalize();
            NCa.normalize();
            CO.normalize();
            OC.normalize();
            CaC.normalize();
            CN1.normalize();
            OCa0.normalize();
            Ca0Ca1.normalize();
        
            nNCaC.normalize();
            nCaCO.normalize();
            nCaCN.normalize();
            nCNCa.normalize();
            nOCN.normalize();            
            nOCa0Ca1.normalize();
            
            dummy = Math.acos(AlgebraFunctions.InnerProduct(nNCaC, nCaCO));

            
            aCaCN = AlgebraFunctions.CalculateAngle(CaC.getInverse(), CN1, nCaCN)*180/Math.PI;
            aCaCO = AlgebraFunctions.CalculateAngle(CaC.getInverse(), CO,  nCaCO)*180/Math.PI;
            aOCN  = AlgebraFunctions.CalculateAngle(OC.getInverse(), CN1, nOCN)*180/Math.PI;
            
            aCNCa = AlgebraFunctions.CalculateAngle(CN0.getInverse(), NCa, nCNCa)*180/Math.PI;
            aNCaC = AlgebraFunctions.CalculateAngle(NCa.getInverse(), CaC, nNCaC)*180/Math.PI;
            aOCa0Ca1 = AlgebraFunctions.CalculateAngle(OCa0.getInverse(), Ca0Ca1, nOCa0Ca1)*180/Math.PI;            
            
            outfile.write((r-1)+","+r+","+aCaCN+","+aCNCa+","+aNCaC+","+aCaCO+","+aOCN+ "," + aOCa0Ca1 + "," + distCO +"\n");
        }
        outfile.close();
    }
    
    public void generateReticle() throws IllegalAccessException
    {
        double maxRadius = 0;
        
        this.lowerBound = new Point3d(Constants.infPos, Constants.infPos, Constants.infPos);
        this.upperBound = new Point3d(Constants.infNeg, Constants.infNeg, Constants.infNeg);

        
        //calculate center and radius for residues, and asigned reticle member
        for(int r=0; r < this.residues.size(); r++)
        {
            this.residues.get(r).calculateCenterRadiusAndBounds();
                    
            if( maxRadius < this.residues.get(r).radius)
            {
                maxRadius = this.residues.get(r).radius;
            }
            
            if( lowerBound.x > this.residues.get(r).lowerBound.x)
            {
                lowerBound.x = this.residues.get(r).lowerBound.x;
            }
            if( upperBound.x < this.residues.get(r).upperBound.x)
            {
                upperBound.x = this.residues.get(r).upperBound.x;
            }

            if( lowerBound.y > this.residues.get(r).lowerBound.y)
            {
                lowerBound.y = this.residues.get(r).lowerBound.y;
            }
            if( upperBound.y < this.residues.get(r).upperBound.y)
            {
                upperBound.y = this.residues.get(r).upperBound.y;
            }

            if( lowerBound.z > this.residues.get(r).lowerBound.z)
            {
                lowerBound.z = this.residues.get(r).lowerBound.z;
            }
            if( upperBound.z < this.residues.get(r).upperBound.z)
            {
                upperBound.z = this.residues.get(r).upperBound.z;
            }
        }
        
        this.reticle = new Reticle(lowerBound.x, upperBound.x, lowerBound.y, upperBound.y, lowerBound.z, upperBound.z);
        this.reticle.maxRadius = maxRadius;
    }
    
    //assign every residue an element in the reticle
    public void assignReticleToResidues(double maxDistance)
    {

        this.reticle.Initialize(maxDistance);
        for( int r = 0 ; r< residues.size(); r++)
        {
        //assign cube
            residues.get(r).reticleCube = reticle.HashCube(residues.get(r).coords.x, residues.get(r).coords.y, residues.get(r).coords.z);
        //add to cube content list
            //initialize cube if it was empty
            if(reticle.cubeMembers.get( residues.get(r).reticleCube ) == null)
            {
                reticle.cubeMembers.set(residues.get(r).reticleCube, new ArrayList<>() );
            }

            reticle.cubeMembers.get( residues.get(r).reticleCube).add(r);
        }
    }

    public void adjustRadius(double factor)
    {
        for (Residue residue : residues) {
            residue.radius *= factor;
        }
    }
    
    public void calculateDihedralAnglesFromCartesian()
    {
        angles.clear();
        Point3d a;
        Point3d b,c,d,e;
        Point3d nAB;
        Point3d nBC, nCD, nDE;
        
        Angles theta;
    //for residue 0:

        theta = new Angles();

    // b = Ca_i    -  N_i
        b = AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.Ca.coords ,  residues.get(0).atomsBackbone.N.coords );
    // c = C_i     -  Ca_i
        c = AlgebraFunctions.SubstractPoints( residues.get(0).atomsBackbone.C.coords ,  residues.get(0).atomsBackbone.Ca.coords );
    // d = N_{i+1} -  C_i 
        d = AlgebraFunctions.SubstractPoints( residues.get(1).atomsBackbone.N.coords ,  residues.get(0).atomsBackbone.C.coords );
    // e = Ca_{i+1}-  N_{i+1}
        e = AlgebraFunctions.SubstractPoints( residues.get(1).atomsBackbone.Ca.coords ,  residues.get(1).atomsBackbone.N.coords );
        
        
//get ortonormal vectors to the planes
        nBC = AlgebraFunctions.CrossProduct(b,c);
        nCD = AlgebraFunctions.CrossProduct(c,d);
        nDE = AlgebraFunctions.CrossProduct(d,e);

        nBC.normalize();
        nCD.normalize();
        nDE.normalize();
        c.normalize();
        d.normalize();
        
        // angle =  acos( - <u,v) / |u|*|v| )
        theta.phi   = 0;
//        theta.psi   = Math.acos( - Functions.InnerProduct(nBC, nCD) );            
        theta.psi   = AlgebraFunctions.CalculateAngle(nBC, nCD, c);
        //theta.omega = Math.acos( - Functions.InnerProduct(nCD, nDE) );
        theta.omega = AlgebraFunctions.CalculateAngle(nCD, nDE, d );

        angles.add(theta);
        
        boolean flagFlip = false;
        
        for(int i = 1; i< residues.size()-1; i++)
        {
            flagFlip = !flagFlip;
            
            theta = new Angles();
            
        // N = 0, Ca = 1, C= 2, O = 3
            
        // a = C_{i-1} - N_i
            a = AlgebraFunctions.SubstractPoints( residues.get(i-1).atomsBackbone.C.coords ,  residues.get(i).atomsBackbone.N.coords );
        // b = N_i     - Ca_i
            b = AlgebraFunctions.SubstractPoints( residues.get( i ).atomsBackbone.N.coords ,  residues.get(i).atomsBackbone.Ca.coords );
        // c = Ca_i    -  C_i
            c = AlgebraFunctions.SubstractPoints( residues.get( i ).atomsBackbone.Ca.coords ,  residues.get(i).atomsBackbone.C.coords );
        // d = C_i     -  N_{i+1} 
            d = AlgebraFunctions.SubstractPoints( residues.get( i ).atomsBackbone.C.coords ,  residues.get(i+1).atomsBackbone.N.coords );
        // e = N_{i+1} -  Ca_{i+1}
            e = AlgebraFunctions.SubstractPoints( residues.get(i+1).atomsBackbone.N.coords ,  residues.get(i+1).atomsBackbone.Ca.coords );
                        
//get ortonormal vectors to the planes

            nAB = AlgebraFunctions.CrossProduct(a,b);
            nBC = AlgebraFunctions.CrossProduct(b,c);
            nCD = AlgebraFunctions.CrossProduct(c,d);
            nDE = AlgebraFunctions.CrossProduct(d,e);

            b.multiplyByScalar(-1);
            c.multiplyByScalar(-1);
            d.multiplyByScalar(-1);

            b.normalize();
            c.normalize();
            d.normalize();

            nAB.normalize();
            nBC.normalize();
            nCD.normalize();
            nDE.normalize();
            
            // angle =  acos( - <u,v) / |u|*|v| )
//            theta.phi   = Math.acos( - Functions.InnerProduct(nAB, nBC) );
//            theta.psi   = Math.acos( - Functions.InnerProduct(nBC, nCD) );            
//            theta.omega = Math.acos( - Functions.InnerProduct(nCD, nDE) );
            theta.phi   = AlgebraFunctions.CalculateAngle(nAB, nBC, b);
            theta.psi   = AlgebraFunctions.CalculateAngle(nBC, nCD, c);            
            theta.omega = AlgebraFunctions.CalculateAngle(nCD, nDE, d);

            angles.add(theta);
        }
        
        theta = new Angles();
        int N = residues.size()-1;
            
        // N = 0, Ca = 1, C= 2, O = 3
            //validate residue positions
            
        // a = C_{i-1} - N_i
            a = AlgebraFunctions.SubstractPoints( residues.get(N-1).atomsBackbone.C.coords ,  residues.get(N).atomsBackbone.N.coords );
        // b = N_i     - Ca_i
            b = AlgebraFunctions.SubstractPoints( residues.get( N ).atomsBackbone.N.coords ,  residues.get(N).atomsBackbone.Ca.coords );
        // c = Ca_i    -  C_i
            c = AlgebraFunctions.SubstractPoints( residues.get( N ).atomsBackbone.Ca.coords ,  residues.get(N).atomsBackbone.C.coords );
                        
//get ortonormal vectors to the planes
            nAB = AlgebraFunctions.CrossProduct(a,b);
            nBC = AlgebraFunctions.CrossProduct(b,c);

            nAB.normalize();
            nBC.normalize();
        
            // angle =  acos( - <u,v) / |u|*|v| )
            //theta.phi   = Math.acos( - Functions.InnerProduct(nAB, nBC) );
            theta.phi   = AlgebraFunctions.CalculateAngle(nAB, nBC, b);
            theta.psi   = 0;            
            theta.omega = 0;

            angles.add(theta);        
    }
    
    public String getDihedralAngles() throws IOException 
    {
        String s = "";
        if(angles.size() ==0)
        {
            calculateCartesianFromDihedral();
        }
        
        for(Angles a : angles)
        {
            s += String.format("%.3f", a.phi) + ","+String.format("%.3f", a.psi) +","+String.format("%.3f", a.omega)+",";
        }
        
        return s.substring(0, s.length()-1);
    }
    
    public void printDihedralAngles(String fileName) throws IOException
    {
        
        BufferedWriter outfile = new BufferedWriter(new FileWriter(fileName));
        outfile.write("id,aminoacid,phi,psi,omega,phi°,psi°,omega°\n");
        for(int i=0; i<this.residues.size()-1;i++)
        {
            outfile.write(i+","+this.residues.get(i).aminoacid+","+angles.get(i).phi+","+angles.get(i).psi+","+angles.get(i).omega);
            outfile.write(","+ (angles.get(i).phi *180/Math.PI) +","+ (angles.get(i).psi *180/Math.PI)+","+(angles.get(i).omega *180/Math.PI));
            outfile.write("\n");
            
        }
        outfile.close();
    }
    
    public void printCoordinatesBackbone(String fileName) throws IOException
    {
        
        BufferedWriter outfile = new BufferedWriter(new FileWriter(fileName));
        outfile.write("id,aminoacid,atom,x,y,z\n");
        for(int i=0; i<this.residues.size()-1;i++)
        {
            outfile.write(this.residues.get(i).id + "," + this.residues.get(i).aminoacid.name + "," + this.residues.get(i).atomsBackbone.N.name  + "," + this.residues.get(i).atomsBackbone.N.coords.x  + "," + this.residues.get(i).atomsBackbone.N.coords.y  + "," + this.residues.get(i).atomsBackbone.N.coords.z+"\n" );
            outfile.write(this.residues.get(i).id + "," + this.residues.get(i).aminoacid.name + "," + this.residues.get(i).atomsBackbone.Ca.name + "," + this.residues.get(i).atomsBackbone.Ca.coords.x + "," + this.residues.get(i).atomsBackbone.Ca.coords.y + "," + this.residues.get(i).atomsBackbone.Ca.coords.z+"\n");
            outfile.write(this.residues.get(i).id + "," + this.residues.get(i).aminoacid.name + "," + this.residues.get(i).atomsBackbone.C.name  + "," + this.residues.get(i).atomsBackbone.C.coords.x  + "," + this.residues.get(i).atomsBackbone.C.coords.y  + "," + this.residues.get(i).atomsBackbone.C.coords.z+"\n" );            
            outfile.write(this.residues.get(i).id + "," + this.residues.get(i).aminoacid.name + "," + this.residues.get(i).atomsBackbone.O.name  + "," + this.residues.get(i).atomsBackbone.O.coords.x  + "," + this.residues.get(i).atomsBackbone.O.coords.y  + "," + this.residues.get(i).atomsBackbone.O.coords.z+"\n" );                        
        }
        outfile.close();
    }
            
    public void test() throws IOException
    {
        Residue r1 = new Residue(1, 'G');
        Residue r2 = new Residue(2, 'A');

        r1.atomsBackbone.N = new Atom(1,"N");
        r1.atomsBackbone.Ca = new Atom(2,"CA");
        r1.atomsBackbone.C = new Atom(3,"C");
        
        r1.atomsBackbone.N.coords.setCoords(0, 0, 0);
        r1.atomsBackbone.N.name = "N";
        
        r1.atomsBackbone.Ca.coords.setCoords(Constants.distNCa, 0, 0);
        r1.atomsBackbone.Ca.name = "Ca";
        
//        r1.atomsBackbone.C.coords.setCoords( Constants.cosAngleNCaC, Constants.sinAngleNCaC, 0);
        
        r1.atomsBackbone.C.coords.setCoords( 1, 0, 0);
        r1.atomsBackbone.C.coords.rotateWithRodriguez(Constants.angleNCaC, new Point3d(0,0,-1));
        
        r1.atomsBackbone.C.coords.multiplyByScalar(Constants.distCaC);
        r1.atomsBackbone.C.coords.sumVector(r1.atomsBackbone.Ca.coords);
        r1.atomsBackbone.C.name = "C";
         
        r1.dihedralAngles.phi = 0;
        r1.dihedralAngles.psi = 3*Math.PI/2;
        r1.dihedralAngles.omega = 0;
        
        r1.vectorNCa = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.Ca.coords, r1.atomsBackbone.N.coords);
        r1.vectorNCa.normalize();
        r1.vectorCaC = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.C.coords, r1.atomsBackbone.Ca.coords);
        r1.vectorCaC.normalize();
        
        //r1.vectorCaC.multiplyByScalar(-1);
        //r1.vectorNCa.multiplyByScalar(-1);
        //r1.calculateNormalVector();

        this.residues.add(r1);

        Angles theta = new Angles();
        boolean flipFlag = true;
/*        theta.omega = 0;
        theta.psi = 0;
        theta.phi = 0;
        r2.setResidue(r1, theta, flipFlag);
        this.residues.add(r2);
*/        
        Residue prev = r1;

        for(int i = 0; i< 4; i++)
        {
            if(flipFlag)
            {
                flipFlag = false;
            }
            else
            {
                flipFlag = true;
            }
            theta.omega = Math.PI/2;
            theta.phi = Math.PI/4;
            theta.psi = - Math.PI/4;
            Residue r = new Residue(i+2, 'A');
            r.setResidue(prev, theta, flipFlag);
            
            prev = r;
            
            this.residues.add(r);
        }
        this.calculateDihedralAnglesFromCartesian();
        this.printDihedralAngles("anglesTest.csv");
/*        
        N2.coords.setCoords(Constants.cosAngleNCCa, Constants.sinAngleNCCa, 0);
        N2.coords.multiplyByScalar(Constants.distNC);
        N2.coords.rotateWithRodriguez(Math.PI/2, Functions.SubstractPoints(C1.coords, Ca1.coords));
        N2.coords.sumVector(C1.coords);
        N2.name = "N";
//        N2.coords.rotateWithRodriguez(0, Functions.SubstractPoints(C1.coords, Ca1.coords));
        
        r2.atoms.add(N2);
        r2.atoms.add(N2);
        r2.atoms.add(N2);
        r2.atoms.add(N2);
        
        this.residues.add(r2);
*/        
        
    }
    
    public void calculateCartesianFromDihedral() throws IOException
    {
        Residue r1;
        if(this.residues.size() != 0)
        {
            r1 = new Residue(1, this.residues.get(0).aminoacid.symbol);
        }
        else
        {
            r1 = new Residue(1, 'A');
        }
          /*  
ATOM      1  N   LEU A  47      64.649  40.034  28.982  1.00 48.72           N  
ATOM      2  CA  LEU A  47      64.804  39.440  30.309  1.00 48.01           C  
ATOM      3  C   LEU A  47      66.268  39.599  30.822  1.00 50.40           C  
ATOM      4  O   LEU A  47      67.200  39.004  30.260  1.00 49.78           O  
    */


        r1.atomsBackbone.N = new Atom(1,"N");
        r1.atomsBackbone.Ca = new Atom(2,"CA");
        r1.atomsBackbone.C = new Atom(3,"C");
        r1.atomsBackbone.O = new Atom(4,"O");
        
        //r1.atomsBackbone.N.coords.setCoords(64.649,  40.034,  28.982);
        r1.atomsBackbone.N.coords.setCoords(31.242,  3.064,  39.284);
        //r1.atomsBackbone.N.name = "N";
        
        //r1.atomsBackbone.Ca.coords.setCoords(64.804,  39.440,  30.309);
        r1.atomsBackbone.Ca.coords.setCoords(31.195,  2.392,  37.963);
        //r1.atomsBackbone.Ca.name = "CA";
        
        //r1.atomsBackbone.C.coords.setCoords(66.268,  39.599,  30.822);
        r1.atomsBackbone.C.coords.setCoords(29.975,  2.923,  37.197);
        //r1.atomsBackbone.C.name = "C";
        
        r1.atomsBackbone.O.coords.setCoords(29.727, 4.132, 37.181);
        //r1.atomsBackbone.O.name = "O";
         
        r1.dihedralAngles.phi = this.angles.get(0).phi;//0;
        r1.dihedralAngles.psi = this.angles.get(0).psi;// 1.9597455914;
        r1.dihedralAngles.omega = this.angles.get(0).omega;//-3.1291640441;

        
        r1.vectorNCa = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.Ca.coords, r1.atomsBackbone.N.coords);
        r1.vectorNCa.normalize();
        r1.vectorCaC = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.C.coords, r1.atomsBackbone.Ca.coords);
        r1.vectorCaC.normalize();
        
        ArrayList<Residue> newResidues = new ArrayList<>();
        newResidues.add(r1);

        Angles theta = new Angles();
        boolean flipFlag = false;
        Residue prev = r1;

        for(int r = 1; r< this.angles.size(); r++)
        {
            flipFlag = !flipFlag;
            theta = this.angles.get(r);
            /*
            theta.omega = q.angles.get(r).omega;
            theta.phi = q.angles.get(r).phi;
            theta.psi = q.angles.get(r).psi;
            */
            Residue res;
            if(this.residues.size() != 0)
            {
                res = new Residue(r+1, this.residues.get(r).aminoacid.name);
            }
            else
            {
                res = new Residue(r+1, "A");
            }
            res.setResidue(prev, theta, flipFlag);
            
            prev = res;
            
            newResidues.add(res);
        }
        //this.calculateDihedralAngles();
        //this.printDihedralAngles(nameFile);
        this.residues.clear();
        this.residues.addAll(newResidues);        
    }

    public void calculateCartesianFromDihedral(String sequence) throws IOException
    {
        Residue r1 = new Residue(1, sequence.charAt(0));
        
          /*  
ATOM      1  N   LEU A  47      64.649  40.034  28.982  1.00 48.72           N  
ATOM      2  CA  LEU A  47      64.804  39.440  30.309  1.00 48.01           C  
ATOM      3  C   LEU A  47      66.268  39.599  30.822  1.00 50.40           C  
ATOM      4  O   LEU A  47      67.200  39.004  30.260  1.00 49.78           O  
    */


        r1.atomsBackbone.N = new Atom(1,"N");
        r1.atomsBackbone.Ca = new Atom(2,"CA");
        r1.atomsBackbone.C = new Atom(3,"C");
        r1.atomsBackbone.O = new Atom(4,"O");
        
        //r1.atomsBackbone.N.coords.setCoords(64.649,  40.034,  28.982);
        r1.atomsBackbone.N.coords.setCoords(31.242,  3.064,  39.284);
        //r1.atomsBackbone.N.name = "N";
        
        //r1.atomsBackbone.Ca.coords.setCoords(64.804,  39.440,  30.309);
        r1.atomsBackbone.Ca.coords.setCoords(31.195,  2.392,  37.963);
        //r1.atomsBackbone.Ca.name = "CA";
        
        //r1.atomsBackbone.C.coords.setCoords(66.268,  39.599,  30.822);
        r1.atomsBackbone.C.coords.setCoords(29.975,  2.923,  37.197);
        //r1.atomsBackbone.C.name = "C";
        
        r1.atomsBackbone.O.coords.setCoords(29.727, 4.132, 37.181);
        //r1.atomsBackbone.O.name = "O";
         
        r1.dihedralAngles.phi = this.angles.get(0).phi;//0;
        r1.dihedralAngles.psi = this.angles.get(0).psi;// 1.9597455914;
        r1.dihedralAngles.omega = this.angles.get(0).omega;//-3.1291640441;

        
        r1.vectorNCa = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.Ca.coords, r1.atomsBackbone.N.coords);
        r1.vectorNCa.normalize();
        r1.vectorCaC = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.C.coords, r1.atomsBackbone.Ca.coords);
        r1.vectorCaC.normalize();
        
        ArrayList<Residue> newResidues = new ArrayList<>();
        newResidues.add(r1);

        Angles theta = new Angles();
        boolean flipFlag = false;
        Residue prev = r1;

        for(int r = 1; r< this.angles.size(); r++)
        {
            flipFlag = !flipFlag;
            theta = this.angles.get(r);
            /*
            theta.omega = q.angles.get(r).omega;
            theta.phi = q.angles.get(r).phi;
            theta.psi = q.angles.get(r).psi;
            */
            Residue res = new Residue(r+1, sequence.charAt(r));
            res.setResidue(prev, theta, flipFlag);
            
            prev = res;
            
            newResidues.add(res);
        }
        //this.calculateDihedralAngles();
        //this.printDihedralAngles(nameFile);
        this.residues.clear();
        this.residues.addAll(newResidues);        
    }

    public void test2(Peptide q, String nameFile) throws IOException
    {
        Residue r1 = new Residue(1, 'L');
        
          /*  
ATOM      1  N   LEU A  47      64.649  40.034  28.982  1.00 48.72           N  
ATOM      2  CA  LEU A  47      64.804  39.440  30.309  1.00 48.01           C  
ATOM      3  C   LEU A  47      66.268  39.599  30.822  1.00 50.40           C  
ATOM      4  O   LEU A  47      67.200  39.004  30.260  1.00 49.78           O  
    */


        r1.atomsBackbone.N = new Atom(1,"N");
        r1.atomsBackbone.Ca = new Atom(2,"CA");
        r1.atomsBackbone.C = new Atom(3,"C");
        
        r1.atomsBackbone.N.coords.setCoords(64.649,  40.034,  28.982);
        r1.atomsBackbone.N.name = "N";
        
        r1.atomsBackbone.Ca.coords.setCoords(64.804,  39.440,  30.309);
        r1.atomsBackbone.Ca.name = "CA";
        
        r1.atomsBackbone.C.coords.setCoords(66.268,  39.599,  30.822);
        r1.atomsBackbone.C.name = "C";
         
        r1.dihedralAngles.phi = q.angles.get(0).phi;//0;
        r1.dihedralAngles.psi = q.angles.get(0).psi;// 1.9597455914;
        r1.dihedralAngles.omega = q.angles.get(0).omega;//-3.1291640441;

        
        r1.vectorNCa = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.Ca.coords, r1.atomsBackbone.N.coords);
        r1.vectorNCa.normalize();
        r1.vectorCaC = AlgebraFunctions.SubstractPoints(r1.atomsBackbone.C.coords, r1.atomsBackbone.Ca.coords);
        r1.vectorCaC.normalize();
        

        this.residues.add(r1);

        Angles theta = new Angles();
        boolean flipFlag = false;
        Residue prev = r1;

        for(int r = 1; r< q.residues.size(); r++)
        {
            flipFlag = !flipFlag;
            theta = q.angles.get(r);
            /*
            theta.omega = q.angles.get(r).omega;
            theta.phi = q.angles.get(r).phi;
            theta.psi = q.angles.get(r).psi;
            */
            Residue res = new Residue(r, q.residues.get(r).aminoacid.name);
            res.setResidue(prev, theta, flipFlag);
            
            prev = res;
            
            this.residues.add(res);
        }
        this.calculateDihedralAnglesFromCartesian();
        this.printDihedralAngles(nameFile);
        
    }
    
    public String getSequence()
    {
        String s = "";
        for(Residue r : residues)
        {
            s += Constants.aminoacidSymbol.get(r.aminoacid.name);
        }
        return s;
    }
    
    public boolean checkFullBackboneReported ()
    {
        for(Residue r : this.residues)
        {
            if(r.atomsBackbone.C == null || r.atomsBackbone.Ca == null || r.atomsBackbone.N == null || r.atomsBackbone.O == null)
            {
                return false;
            }
        }
        return true;
    }
    
    public double RMSDCarbonAlpha(Peptide p)
    {
        double d = 0;
        
        //structures have different size
        if( p.residues.size() != this.residues.size())
        {
            return -1;
        }
        ArrayList<Residue> r = this.residues;
        ArrayList<Residue> s = p.residues;
        
        for( int i=0; i<r.size(); i++)
        {
            d += Math.pow(r.get(i).atomsBackbone.Ca.coords.x - s.get(i).atomsBackbone.Ca.coords.x,2) 
               + Math.pow(r.get(i).atomsBackbone.Ca.coords.y - s.get(i).atomsBackbone.Ca.coords.y,2) 
               + Math.pow(r.get(i).atomsBackbone.Ca.coords.z - s.get(i).atomsBackbone.Ca.coords.z,2);
        }
        return Math.sqrt(d/r.size());
    }
    
    
    public void PrintPDB(String fileName) throws IOException
    {
        BufferedWriter out = new BufferedWriter(new FileWriter(fileName));
        
        String chain = this.chain;
        if(chain==null)
        {
            chain="A";
        }
        
        int atomId=1;
        for(Residue r: this.residues)
        {
//         1         2         3         4         5         6         7         8        
//12345678901234567890123456789012345678901234567890123456789012345678901234567890
//ATOM      1  N   MET A   1      44.061  -3.277   8.755  1.00 22.03           N              
//String.format("%1$10s", title);
            out.write("ATOM  "+                     //1-6   (6) atom
                    String.format("%1$5s",atomId)+  //7-11  (5) atom serial number
                    " " +                           //12    (1) empty space
                    String.format("%1$-4s","N") +    //13-16 (4) atom name
                    " "+                            //17    (1) alternative location indicator
                    String.format("%1$3s", "MET") + //18-20 (3) residue name
                    " " +                           //21    (1) empty space
                    this.chain +                    //22    (1) chain identifier
                    String.format("%1$4s", r.id) +  //23-26 (4) residue sequence number
                    " "  +                          //27    (1) code for insertion of residues
                    "   "+                          //28-30 (3) empty space
                    String.format("%1$8.3f", r.atomsBackbone.N.coords.x) + //31-38    (8) coord x
                    String.format("%1$8.3f", r.atomsBackbone.N.coords.y) + //39-46    (8) coord y
                    String.format("%1$8.3f", r.atomsBackbone.N.coords.z) + //47-54    (8) coord z
                    String.format("%1$6s", "1.00")+ //55-60 (6) occupancy
                    String.format("%1$6s", "22.00")+//61-66 (6) temperature factor
                    "          " +                  //67-76 (10( empty space
                    String.format("%1$2s", "N")+    //77-78 (2) element symbol
                    "  \n");                          //79-80 (2) atom charge
            atomId++;

            out.write("ATOM  "+                     //1-6   (6) atom
                    String.format("%1$5s",atomId)+  //7-11  (5) atom serial number
                    " " +                           //12    (1) empty space
                    String.format("%1$-4s","CA") +    //13-16 (4) atom name
                    " "+                            //17    (1) alternative location indicator
                    String.format("%1$3s", "MET") + //18-20 (3) residue name
                    " " +                           //21    (1) empty space
                    this.chain +                    //22    (1) chain identifier
                    String.format("%1$4s", r.id) +  //23-26 (4) residue sequence number
                    " "  +                          //27    (1) code for insertion of residues
                    "   "+                          //28-30 (3) empty space
                    String.format("%1$8.3f", r.atomsBackbone.Ca.coords.x) + //31-38    (8) coord x
                    String.format("%1$8.3f", r.atomsBackbone.Ca.coords.y) + //39-46    (8) coord y
                    String.format("%1$8.3f", r.atomsBackbone.Ca.coords.z) + //47-54    (8) coord z
                    String.format("%1$6s", "1.00")+ //55-60 (6) occupancy
                    String.format("%1$6s", "22.00")+//61-66 (6) temperature factor
                    "          " +                  //67-76 (10( empty space
                    String.format("%1$2s", "C")+    //77-78 (2) element symbol
                    "  \n");                          //79-80 (2) atom charge
            atomId++;

            out.write("ATOM  "+                     //1-6   (6) atom
                    String.format("%1$5s",atomId)+  //7-11  (5) atom serial number
                    " " +                           //12    (1) empty space
                    String.format("%1$-4s","C") +    //13-16 (4) atom name
                    " " +                           //17    (1) alternative location indicator
                    String.format("%1$3s", "MET") + //18-20 (3) residue name
                    " " +                           //21    (1) empty space
                    this.chain +                    //22    (1) chain identifier
                    String.format("%1$4s", r.id) +  //23-26 (4) residue sequence number
                    " "  +                          //27    (1) code for insertion of residues
                    "   "+                          //28-30 (3) empty space
                    String.format("%1$8.3f", r.atomsBackbone.C.coords.x) + //31-38    (8) coord x
                    String.format("%1$8.3f", r.atomsBackbone.C.coords.y) + //39-46    (8) coord y
                    String.format("%1$8.3f", r.atomsBackbone.C.coords.z) + //47-54    (8) coord z
                    String.format("%1$6s", "1.00")+ //55-60 (6) occupancy
                    String.format("%1$6s", "22.00")+//61-66 (6) temperature factor
                    "          " +                  //67-76 (10( empty space
                    String.format("%1$2s", "C")+    //77-78 (2) element symbol
                    "  \n");                          //79-80 (2) atom charge
            atomId++;

            out.write("ATOM  "+                     //1-6   (6) atom
                    String.format("%1$5s",atomId)+  //7-11  (5) atom serial number
                    " " +                           //12    (1) empty space
                    String.format("%1$-4s","O") +    //13-16 (4) atom name
                    " "+                            //17    (1) alternative location indicator
                    String.format("%1$3s", "MET") + //18-20 (3) residue name
                    " " +                           //21    (1) empty space
                    this.chain +                    //22    (1) chain identifier
                    String.format("%1$4s", r.id) +  //23-26 (4) residue sequence number
                    " "  +                          //27    (1) code for insertion of residues
                    "   "+                          //28-30 (3) empty space
                    String.format("%1$8.3f", r.atomsBackbone.O.coords.x) + //31-38    (8) coord x
                    String.format("%1$8.3f", r.atomsBackbone.O.coords.y) + //39-46    (8) coord y
                    String.format("%1$8.3f", r.atomsBackbone.O.coords.z) + //47-54    (8) coord z
                    String.format("%1$6s", "1.00")+ //55-60 (6) occupancy
                    String.format("%1$6s", "22.00")+//61-66 (6) temperature factor
                    "          " +                  //67-76 (10( empty space
                    String.format("%1$2s", "O")+    //77-78 (2) element symbol
                    "  \n");                          //79-80 (2) atom charge
            atomId++;
        }
        out.close();
    }
    
    public void loadCif (String fileName, String chain, boolean flagBackbone) throws FileNotFoundException, IOException
    {
        BufferedReader in;
        
        if(fileName.contains(".gz"))
        {
            //GZIPInputStream gzip = new GZIPInputStream(new FileInputStream("F:/gawiki-20090614-stub-meta-history.xml.gz"));
            in = new BufferedReader(new InputStreamReader( new GZIPInputStream(new FileInputStream(fileName)) ));
            //infile = new BufferedReader(new Gz FileReader(fileName));            
        }
        else
        {
            in = new BufferedReader(new FileReader(fileName));
        }
        
        double x,y,z;
        this.chain = chain;
        String line = "";
        String atomName = "";
        String resName = "";
        String chainID = "";
        Integer resSeq = 0;
        Integer serial = 0;
        
        while( (line = in.readLine())!= null)
        {
            if(! line.startsWith("ATOM "))
            {
                continue;
            }

            while(line.contains("  "))
            {
                line = line.replaceAll("  ", " ");
            }
            String[] fields = line.split(" ");
            
            /*
        0     _atom_site.group_PDB 
        1     _atom_site.id 
        2     _atom_site.type_symbol 
        3     _atom_site.label_atom_id 
        4     _atom_site.label_alt_id 
        5     _atom_site.label_comp_id 
        6     _atom_site.label_asym_id 
        7     _atom_site.label_entity_id 
        8     _atom_site.label_seq_id 
        9     _atom_site.pdbx_PDB_ins_code 
        10    _atom_site.Cartn_x 
        11    _atom_site.Cartn_y 
        12    _atom_site.Cartn_z 
        13    _atom_site.occupancy 
        14    _atom_site.B_iso_or_equiv 
        15    _atom_site.pdbx_formal_charge 
        16    _atom_site.auth_seq_id 
        17    _atom_site.auth_comp_id 
        18    _atom_site.auth_asym_id 
        19    _atom_site.auth_atom_id 
        20    _atom_site.pdbx_PDB_model_num 
            0      1      2  3     4 5   6  7  8   9 10      11      12      13   14    15 16  17  18 19    20
            ATOM   86819  N  N     . THR OB 66 1   ? 240.774 299.008 273.887 1.00 85.41  ? 2   THR AJ N     1 
            ATOM   86820  C  CA    . THR OB 66 1   ? 239.733 299.148 272.880 1.00 82.21  ? 2   THR AJ CA    1 
            ATOM   86821  C  C     . THR OB 66 1   ? 240.312 299.572 271.536 1.00 79.25  ? 2   THR AJ C     1 
            ATOM   86822  O  O     . THR OB 66 1   ? 240.452 300.765 271.268 1.00 79.09  ? 2   THR AJ O     1 
            ATOM   86823  C  CB    . THR OB 66 1   ? 238.670 300.167 273.310 1.00 20.00  ? 2   THR AJ CB    1 

            */
            atomName = fields[3];
            
            if( !flagBackbone )
            {
                // skip non-backbone atoms
                if( !atomName.equals("C") && !atomName.equals("O") && !atomName.equals("N") && !atomName.equals("CA"))
                {
                    continue;
                }
            }
            
                //aminoacid name: PRO, GLN, VAL, ...
            resName = fields[5];
                //residue number 
            resSeq = Integer.parseInt(fields[8]);
                //atom number
            /*
            if(flagBackbone)
             {
              if(atomName.equals("C") || atomName.equals("O") || atomName.equals("N") || atomName.equals("CA"))
                 serial = Integer.parseInt(line.substring(6, 11).trim());
             }
            else serial = Integer.parseInt(line.substring(6, 11).trim());
            */
            serial = Integer.parseInt(fields[1]);
            x = Double.parseDouble(fields[10]);
            y = Double.parseDouble(fields[11]);
            z = Double.parseDouble(fields[12]);
            if(!serial.equals(""))
             {
                //search can be skipped if the .pdb file is sorted by residue id
                Residue currentResidue;
                int residueId = this.residuesInitialized.indexOf(resSeq);
                if(residueId == -1)
                {
                    this.residuesInitialized.add(resSeq);
                    currentResidue = new Residue(resSeq, resName);
                    //currentResidue.id = resSeq;
                    //currentResidue.aminoacid = resName;
                    //currentResidue.atoms = new ArrayList<>();
                    this.residues.add(currentResidue);
                }
                else
                {
                    currentResidue = this.residues.get(residueId);
                }
                //add atom to residue
                Atom a = new Atom(serial, atomName, x, y, z);
                if(atomName.contentEquals("N") )
                {
                    currentResidue.atomsBackbone.N = a;
                }
                else if(atomName.contentEquals("C") )
                {
                    currentResidue.atomsBackbone.C = a;
                }
                else if(atomName.contentEquals("CA") )
                {
                    currentResidue.atomsBackbone.Ca = a;
                }
                else if(atomName.contentEquals("O"))
                {
                    currentResidue.atomsBackbone.O = a;
                }
                else
                {
                    currentResidue.atomsLateralChain.add(a);
                }
             }            
        }
    }
}
