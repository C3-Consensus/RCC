
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.NoSuchElementException;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author fonty
 */
public class Residue
{
    int id;
    Aminoacid aminoacid;
    // N = 0, Ca = 1, C= 2, O = 3
    Backbone atomsBackbone;
    ArrayList<Atom> atomsLateralChain;
    ArrayList<Integer> atomsIds;
    public Point3d coords;
    public double radius;
    public int reticleCube;
    public Angles dihedralAngles;
    // CaN = Ca - N
    Point3d vectorNCa;
    // CCa = C  - Ca
    Point3d vectorCaC;
    // n = normal vector to the plane N, Ca, C
    Point3d vectorNormal;
    //residue is contained in the cube formed by the upper and lower bounds
    Point3d lowerBound;
    Point3d upperBound;
    
    public Residue()
    {
        this.atomsLateralChain = new ArrayList<>();
        this.atomsIds = new ArrayList<>();
        this.dihedralAngles = new Angles();
        this.vectorCaC = null;
        this.vectorNCa = null;
        this.vectorNormal = null;
        this.coords = new Point3d();
        this.atomsBackbone = new Backbone();       
    }
    
    public Residue(int id, String name)
    {
        this.id = id;
        this.aminoacid = new Aminoacid(name);
        this.atomsLateralChain = new ArrayList<>();
        this.atomsIds = new ArrayList<>();
        this.dihedralAngles = new Angles();
        this.vectorCaC = null;
        this.vectorNCa = null;
        this.vectorNormal = null;
        this.coords = new Point3d();
        this.atomsBackbone = new Backbone();
    }
    
    public Residue(int id, char symbol)
    {
        this.id = id;
        this.aminoacid = new Aminoacid(symbol);
        this.atomsLateralChain = new ArrayList<>();
        this.atomsIds = new ArrayList<>();
        this.dihedralAngles = new Angles();
        this.vectorCaC = null;
        this.vectorNCa = null;
        this.vectorNormal = null;
        this.coords = new Point3d();
        this.atomsBackbone = new Backbone();
    }
    
    public Residue(Residue a)
    {
        id = a.id;
        aminoacid = new Aminoacid(a.aminoacid);
        atomsBackbone = new Backbone(a.atomsBackbone);
        atomsLateralChain = new ArrayList<>();
        for(Atom at : a.atomsLateralChain)
        {
            atomsLateralChain.add(new Atom(at));
        }
        atomsIds = new ArrayList<>();
        for(Integer i : a.atomsIds)
        {
            atomsIds.add(new Integer(i));
        }
        coords = new Point3d(a.coords);
        radius = a.radius;
        reticleCube = a.reticleCube;
        dihedralAngles = new Angles(a.dihedralAngles);
        
        if(a.vectorNCa != null)
        {
            vectorNCa = new Point3d(a.vectorNCa);
        }
        else
        {
            vectorNCa = null;
        }
        if(a.vectorCaC != null)
        {
            vectorCaC = new Point3d(a.vectorCaC);
        }
        else
        {
            vectorCaC = null;
        }
        
        if(a.vectorNormal != null)
        {
            vectorNormal = new Point3d(a.vectorNormal);
        }
        else
        {
            vectorNormal = null;
        }
        
        if(a.lowerBound !=null)
        {
            lowerBound = new Point3d(a.lowerBound);
        }
        else
        {
            lowerBound = null;
        }
        if(a.upperBound != null)
        {
            upperBound = new Point3d(a.upperBound);
        }
        else
        {
            upperBound = null;
        }
    }
    
    public Iterator<Atom> iterator(boolean backboneOnly)
    {
        if(backboneOnly)
        {
            return new residueIteratorBackbone();
        }
        return new residueIterator();
    }
    
    private class residueIterator implements Iterator<Atom> {
        private Atom cursor;
        private int indexBackbone;
        private int indexLateralChain;
        private int totalAtoms;
        private int indexCurrent;
        private int nullAtomsBackbone;

        public residueIterator() {
            this.cursor = null;
            nullAtomsBackbone = 0;
            if( Residue.this.atomsBackbone.N == null)
            {
                nullAtomsBackbone++;
            }
            if( Residue.this.atomsBackbone.C == null)
            {
                nullAtomsBackbone++;
            }
            if( Residue.this.atomsBackbone.Ca == null)
            {
                nullAtomsBackbone++;
            }
            if( Residue.this.atomsBackbone.O == null)
            {
                nullAtomsBackbone++;
            }
            
            totalAtoms = Residue.this.atomsLateralChain.size() + 4 - nullAtomsBackbone;
            indexBackbone = 0;
            indexLateralChain =0;
        }

        public boolean hasNext() {
            return this.indexCurrent < this.totalAtoms;
        }

        public Atom next() {
            if(indexBackbone < 4)
            {
                if(indexBackbone == 0)
                {
                    this.cursor = Residue.this.atomsBackbone.N;
                    indexBackbone++;
                    if(this.cursor != null)
                    {
                        indexCurrent++;
                        return this.cursor;
                    }
                }
                if(indexBackbone == 1)
                {
                    this.cursor = Residue.this.atomsBackbone.Ca;
                    indexBackbone++;
                    if(this.cursor != null)
                    {
                        indexCurrent++;
                        return this.cursor;
                    }
                }
                if(indexBackbone == 2)
                {
                    this.cursor = Residue.this.atomsBackbone.C;
                    indexBackbone++;
                    if(this.cursor != null)
                    {
                        indexCurrent++;
                        return this.cursor;
                    }
                }
                if(indexBackbone == 3)
                {
                    this.cursor = Residue.this.atomsBackbone.O;
                    indexBackbone++;
                    if(this.cursor != null)
                    {
                        indexCurrent++;
                        return this.cursor;
                    }
                }
            }
//            else
//            {
            //if(indexLateralChain < residue.this.atomsLateralChain.size())
                
                this.cursor = Residue.this.atomsLateralChain.get(indexLateralChain);
                indexLateralChain++;
//            }
            indexCurrent++;
            return this.cursor;
            //throw new NoSuchElementException();
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
    
    private class residueIteratorBackbone implements Iterator<Atom> {
        private Atom cursor;
        private int indexBackbone;
        private int totalAtoms;

        public residueIteratorBackbone() {
            this.cursor = null;
            totalAtoms =  0;
            if( Residue.this.atomsBackbone.N != null)
            {
                totalAtoms++;
            }
            if( Residue.this.atomsBackbone.C != null)
            {
                totalAtoms++;
            }
            if( Residue.this.atomsBackbone.Ca != null)
            {
                totalAtoms++;
            }
            if( Residue.this.atomsBackbone.O != null)
            {
                totalAtoms++;
            }
            
            
            indexBackbone = 1;
        }

        public boolean hasNext() {
            return this.indexBackbone < this.totalAtoms;
        }
@Deprecated
        public Atom nextOld() {
            if(indexBackbone < 4)
            {
                if(indexBackbone == 0)
                {
                    this.cursor = Residue.this.atomsBackbone.N;
                    indexBackbone++;
                    return this.cursor;
                }
                if(indexBackbone == 1)
                {
                    this.cursor = Residue.this.atomsBackbone.Ca;
                    indexBackbone++;
                    return this.cursor;
                }
                if(indexBackbone == 2)
                {
                    this.cursor = Residue.this.atomsBackbone.C;
                    indexBackbone++;
                    return this.cursor;
                }
                if(indexBackbone == 3)
                {
                    this.cursor = Residue.this.atomsBackbone.O;
                    indexBackbone++;
                    return this.cursor;
                }
            }
            
            throw new NoSuchElementException();
        }
        public Atom next() {
            if(indexBackbone < 4)
            {
                if(indexBackbone == 0)
                {
                    this.cursor = Residue.this.atomsBackbone.N;
                    indexBackbone++;
                    if(this.cursor!=null)
                    {
                        return this.cursor;
                    }
                }
                if(indexBackbone == 1)
                {
                    this.cursor = Residue.this.atomsBackbone.Ca;
                    indexBackbone++;
                    if(this.cursor!=null)
                    {
                        return this.cursor;
                    }
                }
                if(indexBackbone == 2)
                {
                    this.cursor = Residue.this.atomsBackbone.C;
                    indexBackbone++;
                    if(this.cursor!=null)
                    {
                        return this.cursor;
                    }
                }
                if(indexBackbone == 3)
                {
                    this.cursor = Residue.this.atomsBackbone.O;
                    indexBackbone++;
                    if(this.cursor!=null)
                    {
                        return this.cursor;
                    }
                }
            }
//            return null;
            throw new NoSuchElementException();
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }
    
    
    public void calculateCenterRadiusAndBounds() throws IllegalArgumentException, IllegalAccessException 
    {
        this.lowerBound = new Point3d(Constants.infPos, Constants.infPos, Constants.infPos);
        this.upperBound = new Point3d(Constants.infNeg, Constants.infNeg, Constants.infNeg);
        /*
        double minX = Constants.infPos;
        double maxX = Constants.infNeg;
        double minY = Constants.infPos;
        double maxY = Constants.infNeg;
        double minZ = Constants.infPos;
        double maxZ = Constants.infNeg;
        */
        this.coords.x = 0;
        this.coords.y = 0;
        this.coords.z = 0;
    
    //calculate center
        int numFields = 0;
        for(Field f : this.atomsBackbone.getClass().getFields())
        {
            if(f.get(this.atomsBackbone)==null)
            {
                continue;
            }
            numFields++;
        //sum to calculate center of mass
            coords.sumVector( ((Atom) f.get(this.atomsBackbone) ).coords );
        //calculate bounds
            if( ((Atom) f.get(this.atomsBackbone) ).coords.x > upperBound.x )
            {
                upperBound.x = ((Atom) f.get(this.atomsBackbone) ).coords.x;
            }
            if( ((Atom) f.get(this.atomsBackbone) ).coords.x < lowerBound.x )
            {
                lowerBound.x = ((Atom) f.get(this.atomsBackbone) ).coords.x;                
            }
            
            if( ((Atom) f.get(this.atomsBackbone) ).coords.y > upperBound.y )
            {
                upperBound.y = ((Atom) f.get(this.atomsBackbone) ).coords.y;                
            }
            if( ((Atom) f.get(this.atomsBackbone) ).coords.y < lowerBound.y )
            {
                lowerBound.y = ((Atom) f.get(this.atomsBackbone) ).coords.y;
            }
            
            if( ((Atom) f.get(this.atomsBackbone) ).coords.z > upperBound.z )
            {
                upperBound.z = ((Atom) f.get(this.atomsBackbone) ).coords.z;                
            }
            if( ((Atom) f.get(this.atomsBackbone) ).coords.z < lowerBound.z )
            {
                lowerBound.z = ((Atom) f.get(this.atomsBackbone) ).coords.z;                
            }
            
        }
        
        for(int at = 0; at < this.atomsLateralChain.size(); at++)
        {
            this.coords.x += this.atomsLateralChain.get(at).coords.x;
            this.coords.y += this.atomsLateralChain.get(at).coords.y;
            this.coords.z += this.atomsLateralChain.get(at).coords.z;
            
            lowerBound.x = Math.min(lowerBound.x, this.atomsLateralChain.get(at).coords.x);
            lowerBound.y = Math.min(lowerBound.y, this.atomsLateralChain.get(at).coords.y);
            lowerBound.z = Math.min(lowerBound.z, this.atomsLateralChain.get(at).coords.z);
            
            upperBound.x = Math.max(upperBound.x, this.atomsLateralChain.get(at).coords.x);
            upperBound.y = Math.max(upperBound.y, this.atomsLateralChain.get(at).coords.y);
            upperBound.z = Math.max(upperBound.z, this.atomsLateralChain.get(at).coords.z);
            
        }
        this.coords.x /= numFields + atomsLateralChain.size();
        this.coords.y /= numFields + atomsLateralChain.size();
        this.coords.z /= numFields + atomsLateralChain.size();
        //calculate radius
        this.radius=0;
        
        for(int at=0;at < this.atomsLateralChain.size(); at++)
        {
            double dist =  AlgebraFunctions.calculateDistance(this.coords, this.atomsLateralChain.get(at).coords);
            if( dist > this.radius)
            {
                this.radius = dist;
            }
        }
        for(Field f : this.atomsBackbone.getClass().getFields())
        {
            Atom backboneAtom = (Atom) f.get(this.atomsBackbone);
            double dist =  AlgebraFunctions.calculateDistance(this.coords, backboneAtom.coords);
            if( dist > this.radius)
            {
                this.radius = dist;
            }
        }
    }
    
    public void calculateNormalVector(boolean flipFlag)
    {
        if( this.vectorNormal != null)
        {
            return;
        }
        
        if(flipFlag)
        {
            this.vectorNormal = AlgebraFunctions.CrossProduct(this.vectorCaC, this.vectorNCa);
        }
        else
        {
            this.vectorNormal = AlgebraFunctions.CrossProduct(this.vectorNCa, this.vectorCaC);
        }
        this.vectorNormal.normalize();
    }
    
    
    public void setResidue( Residue anterior, Angles theta, boolean flipFlag  )
    {
        Point3d vectorCN;
        Point3d vectorNormalCaCN, vectorNormalCNCa;
        //Angles must be flipped by 360 every residue
//set N
    // rotate rotation on plane with previous residue according to bond angle
    // rotate with respect to dihedral angle
    // rescale vector
    // move vector to position
        
        anterior.calculateNormalVector(false);
    //set using saved normalized vector from previous residue
        atomsBackbone.N = new Atom(1, "N", anterior.vectorCaC);
    //rotate by its standard angle
        if(flipFlag)
        {
            atomsBackbone.N.coords.rotateWithRodriguez(Constants.angleCaCNFlip, anterior.vectorNormal);
        }
        else
        {
            atomsBackbone.N.coords.rotateWithRodriguez(Constants.angleCaCN, anterior.vectorNormal);            
        }
    //rotate by dihedral angle
        atomsBackbone.N.coords.rotateWithRodriguez(anterior.dihedralAngles.psi, anterior.vectorCaC.getInverse());
        //atomsBackbone.N.coords.rotateWithRodriguez(anterior.dihedralAngles.psi, anterior.vectorCaC);
    //save normalized vector for futur atom position calculations
        vectorCN = new Point3d(atomsBackbone.N.coords);
    //rescalel vector to its standard size
        atomsBackbone.N.coords.multiplyByScalar(Constants.distCN);
    //shift vector to its correct position in the structure
        atomsBackbone.N.coords.sumVector(anterior.atomsBackbone.C.coords);

//set Ca
        atomsBackbone.Ca = new Atom(2, "Ca",vectorCN);
//        vectorNormalCaCN = Functions.CrossProductNormalized(anterior.vectorCaC, vectorCN);
        if(flipFlag)
        {
            vectorNormalCaCN = AlgebraFunctions.CrossProductNormalized(anterior.vectorCaC, vectorCN);
            atomsBackbone.Ca.coords.rotateWithRodriguez( Constants.angleCNCaFlip , vectorNormalCaCN);
        }
        else
        {
            vectorNormalCaCN = AlgebraFunctions.CrossProductNormalized(vectorCN, anterior.vectorCaC);
            atomsBackbone.Ca.coords.rotateWithRodriguez( Constants.angleCNCa , vectorNormalCaCN);
            
        }
        atomsBackbone.Ca.coords.rotateWithRodriguez(anterior.dihedralAngles.omega , vectorCN.getInverse());
        //atomsBackbone.Ca.coords.rotateWithRodriguez(anterior.dihedralAngles.omega , vectorCN);
        this.vectorNCa = new Point3d(atomsBackbone.Ca.coords);        
        atomsBackbone.Ca.coords.multiplyByScalar(Constants.distNCa);
        atomsBackbone.Ca.coords.sumVector( atomsBackbone.N.coords );
        
        
//set C
        atomsBackbone.C = new Atom(3, "C",this.vectorNCa);
        //vectorNormalCNCa = Functions.CrossProductNormalized(vectorCN, this.vectorNCa);
        if(flipFlag)
        {
            vectorNormalCNCa = AlgebraFunctions.CrossProductNormalized(vectorCN, this.vectorNCa);
            atomsBackbone.C.coords.rotateWithRodriguez( Constants.angleNCaCFlip , vectorNormalCNCa);
        }
        else
        {
            vectorNormalCNCa = AlgebraFunctions.CrossProductNormalized(this.vectorNCa, vectorCN);
            atomsBackbone.C.coords.rotateWithRodriguez( Constants.angleNCaC , vectorNormalCNCa);
        }
        atomsBackbone.C.coords.rotateWithRodriguez( theta.phi , this.vectorNCa.getInverse());
        this.vectorCaC = new Point3d(atomsBackbone.C.coords);
        atomsBackbone.C.coords.multiplyByScalar(Constants.distCaC);
        atomsBackbone.C.coords.sumVector(atomsBackbone.Ca.coords);
        
//set O for previous residue
        boolean flagSetOxygenWithNextCA = false;
        if(flagSetOxygenWithNextCA)
        {
            Point3d CaCa = new  Point3d(this.atomsBackbone.Ca.coords, anterior.atomsBackbone.Ca.coords);
            CaCa.normalize();

            Point3d vectorCa0N = new Point3d(atomsBackbone.N.coords, anterior.atomsBackbone.Ca.coords);
            vectorCa0N.normalize();           

            Point3d vectorNormalCa0NCa = AlgebraFunctions.CrossProductNormalized(vectorCa0N, CaCa);

            Point3d[] posibleOxygen = setOxygenNextCa(vectorNormalCa0NCa, Constants.distCaO[anterior.aminoacid.numAminoacid], Constants.angleCa0OCa1Flip[anterior.aminoacid.numAminoacid], CaCa);

            posibleOxygen[0].sumVector(anterior.atomsBackbone.Ca.coords);
            posibleOxygen[1].sumVector(anterior.atomsBackbone.Ca.coords);

            double dist1 = AlgebraFunctions.calculateDistance( posibleOxygen[0], anterior.atomsBackbone.C.coords);
            double dist2 = AlgebraFunctions.calculateDistance( posibleOxygen[1], anterior.atomsBackbone.C.coords);

            //correct position has to be closer to the Carbon atom
            if( dist1 < dist2)
            {
                anterior.atomsBackbone.O = new Atom(4, "O", posibleOxygen[0]);
            }
            else
            {
                anterior.atomsBackbone.O = new Atom(4, "O", posibleOxygen[1]);
            }
        }
        
/***/
  //set vector as Ca - Ca_-1
        else
        {
            anterior.atomsBackbone.O = new Atom(4, "O", anterior.vectorCaC);
            this.calculateNormalVector(false);
/**/
            Point3d CaCa = new  Point3d(anterior.atomsBackbone.C.coords, anterior.atomsBackbone.Ca.coords);
            CaCa.normalize();

            Point3d vectorCa0N = new Point3d(atomsBackbone.N.coords, anterior.atomsBackbone.Ca.coords);
            //            point3d vectorCa0N = new point3d(atomsBackbone.Ca.coords, anterior.atomsBackbone.Ca.coords);
            vectorCa0N.normalize();           

            Point3d vectorNormalCa0NCa = AlgebraFunctions.CrossProductNormalized(vectorCa0N, CaCa);
/**/
        //rotate by its standard angle
            if(!flipFlag)
            {
                anterior.atomsBackbone.O.coords.rotateWithRodriguez(Constants.angleCaCOFlip, vectorNormalCa0NCa);
            }
            else
            {
                anterior.atomsBackbone.O.coords.rotateWithRodriguez(Constants.angleCaC0, vectorNormalCa0NCa.getInverse());            
            }
        // no dihedral angle to rotate by
        // no need to save normalized vector for futur atom position calculations
        //rescalel vector to its standard size
            anterior.atomsBackbone.O.coords.multiplyByScalar(Constants.distCO);
        //shift vector to its correct position in the structure
            anterior.atomsBackbone.O.coords.sumVector(anterior.atomsBackbone.C.coords);            

        }
                
                
        this.dihedralAngles = theta;
        this.vectorCaC.normalize();
        this.vectorNCa.normalize();
        
    }
    
    public static Point3d[] setOxygenNextCa(Point3d normalVector, double dist, double angle, Point3d u )
    {
        double x1,y1,z1,x2,y2,z2;
        double A,B,C,D,E,F,G,H,I;
        
        double tt = normalVector.getNorm();
        double uu = u.getNorm();
        
        A = normalVector.x;
        B = normalVector.y;
        C = normalVector.z;
        /*
        D = (A * u.y) - (B * u.x);
        E = A * dist * Math.cos(angle);
        F = (C * u.x) - (A * u.z);
        
        double DCBF = D*C + B*F;
        double A2 = A * A;
        double A2D2 = A2 * D*D;
        
        G = A2D2 + F*F*A2 + DCBF*DCBF;
        H = 2* ( A2*E*F + DCBF* B *E );
        I = E*E * (A2 + B * B) - (A2D2 * dist * dist);
        
        double disc = Math.sqrt(  H*H - 4 * G * I);
        double GG = 2*G;
        
        z1 = (-H + disc)/GG;
        z2 = (-H - disc)/GG;
        
        x1 = - normalVector.z * z1 - normalVector.y * ( E + F * z1 )/ D;
        x2 = - normalVector.z * z2 - normalVector.y * ( E + F * z2 )/ D;
        
        y1 = - ( normalVector.x * x1 + normalVector.z * z1 ) / normalVector.y;
        y2 = - ( normalVector.x * x2 + normalVector.z * z2 ) / normalVector.y;
        /**/
        
        A = normalVector.x;
        B = normalVector.y;
        C = normalVector.z;
        
        D = A * u.y - B * u.x;
        E = A * Math.cos(angle) * dist;
        F = A * u.z - C * u.x;
        
        G = A*A * D*D + F*F * A*A + Math.pow(B*F - D*C ,2);
        H = - 2*B*E* (B*F - D*C) - 2 * E * F *A*A;
        I = B*B * E*E + E*E * A*A - A*A * D*D * dist*dist;
        
        double disc = Math.sqrt(H*H - 4*G*I);
        
        z1 = (-H + disc ) / (2*G);
        z2 = (-H - disc ) / (2*G);
        
        
        double dum1 = G*z1*z1 + H*z1 + I;
        double dum2 = G*z2*z2 + H*z2 + I;
        
        
        y1 = (E - F*z1 ) / D;
        y2 = (E - F*z2 ) / D;
        x1 = - (B* y1 / A) - (C * z1)/A;
        x2 = - (B* y2 / A) - (C * z2)/A;
        
/*       
        point3d p1 = new point3d(x1,y1,z1);
        point3d p2 = new point3d(x2,y2,z2);
        
        double dummy = Functions.InnerProduct(normalVector, p1);
        double dummy2 = Functions.InnerProduct(normalVector, p2);
        
        double dummy3 = Functions.InnerProduct(u, p1)/p1.getNorm();
        double dummy4 = Functions.InnerProduct(u, p2)/p2.getNorm();
        
        double dummy5 = p1.getNorm();
        double dummy6 = p2.getNorm();
        
        double dummy7 = Math.cos(angle);
*/        
        Point3d[] v = new Point3d[2];
        v[0] = new Point3d(x1, y1, z1);
        v[1] = new Point3d(x2, y2, z2);
        return v;
    }
}
