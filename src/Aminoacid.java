/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Fonty
 */
public class Aminoacid {
    String name;
    char symbol;
    int numAminoacid;
    public Aminoacid(String name)
    {
        this.name = name;
        switch (name.toUpperCase())
        {
            case "ALA":   symbol = 'a';   break;
            case "CYS":   symbol = 'c';   break;
            case "ASP":   symbol = 'd';   break;
            case "GLU":   symbol = 'e';   break;
            case "PHE":   symbol = 'f';   break;
            case "GLY":   symbol = 'g';   break;
            case "HIS":   symbol = 'h';   break;
            case "ILE":   symbol = 'i';   break;
            case "LYS":   symbol = 'k';   break;
            case "LEU":   symbol = 'l';   break;
            case "MET":   symbol = 'm';   break;
            case "ASN":   symbol = 'n';   break;
            case "PRO":   symbol = 'p';   break;
            case "GLN":   symbol = 'q';   break;
            case "ARG":   symbol = 'r';   break;
            case "SER":   symbol = 's';   break;
            case "THR":   symbol = 't';   break;
            case "VAL":   symbol = 'v';   break;
            case "TRP":   symbol = 'w';   break;
            case "TYR":   symbol = 'y';   break;
            default:      symbol = 0;     break;
        }
        numAminoacid = symbol - 'a';
    }
    
    public Aminoacid(char symbol)
    {
        if(symbol < 'a')
        {
            symbol = (char) (symbol -'A' + 'a');
        }
        this.symbol = symbol;
        numAminoacid = symbol - 'a';
        
        switch(symbol)
        {
            case 'a':	this.name = "ALA";	break;
            case 'c':	this.name = "CYS";	break;
            case 'd':	this.name = "ASP";	break;
            case 'e':	this.name = "GLU";	break;
            case 'f':	this.name = "PHE";	break;
            case 'g':	this.name = "GLY";	break;
            case 'h':	this.name = "HIS";	break;
            case 'i':	this.name = "ILE";	break;
            case 'k':	this.name = "LYS";	break;
            case 'l':	this.name = "LEU";	break;
            case 'm':	this.name = "MET";	break;
            case 'n':	this.name = "ASN";	break;
            case 'p':	this.name = "PRO";	break;
            case 'q':	this.name = "GLN";	break;
            case 'r':	this.name = "ARG";	break;
            case 's':	this.name = "SER";	break;
            case 't':	this.name = "THR";	break;
            case 'v':	this.name = "VAL";	break;
            case 'w':	this.name = "TRP";	break;
            case 'y':	this.name = "TYR";	break;
            default:	this.name = "invalid";  break;
            
        }
    }
    
    public Aminoacid(Aminoacid a)
    {
        this.name = new String(a.name);
        this.symbol = a.symbol;
        this.numAminoacid = a.numAminoacid;
    }
}
