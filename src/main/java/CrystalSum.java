import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import java.util.ArrayList;

public class CrystalSum {
    double a=1;
    double c=2.077294686;
    int num_in_cell=4;
    RealVector[] latticeVectors = new RealVector[3];
    double[][] basis;
    {
        latticeVectors[0] = new ArrayRealVector(new double[]{a, 0, 0});
        latticeVectors[1] = new ArrayRealVector(new double[]{0, a, 0});
        latticeVectors[2] = new ArrayRealVector(new double[]{0, 0, c});
         basis = new double[][]{
                {0, 0, 0},
                {0.5, 0, 0.25},
                {0.5, 0.5, 0.5},
                {0, 0.5, 0.75}
        };
    }
    public double dipolar(int comp1, int comp2, RealVector ri, RealVector rj){
        double norm = ri.subtract(rj).getNorm();
        if (comp1==comp2){
            return (norm*norm - 3*Math.pow(ri.subtract(rj).getEntry(comp1), 2))/Math.pow(norm, 5);
        } else {
            return (-3*ri.subtract(rj).getEntry(comp1)*ri.subtract(rj).getEntry(comp2))/Math.pow(norm, 5);
        }
    }

    public RealVector getLocation(int n, int[] limitingBox){
        int Lx = limitingBox[0], Ly = limitingBox[1];
        int k = ((n/num_in_cell)/Lx)/Ly/8;
        int sym_op = ((n/num_in_cell)/Lx)/Ly%8;
        int j = ((n/num_in_cell)/Lx)%Ly;
        int i = (n/num_in_cell)%Lx;
        int l = n%num_in_cell;

        int[] symVec = new int[]{sym_op%2,sym_op/2%2,sym_op/4%2};

        int[] locCoordinates = new int[]{i, j, k};
        RealVector loc = new ArrayRealVector(3);
        for (int coor=0;coor<latticeVectors.length;coor++){
            try {
                loc.combineToSelf(1, locCoordinates[coor], latticeVectors[coor]);   // add primitive lattice vectors according to (i,j,k)
                loc.combineToSelf(1, basis[l][coor], latticeVectors[coor]);         // add fraction of primitive lattice vector according to basis index l
                if (coor != 2 && symVec[coor] != 0) {
                    loc.combineToSelf(1, -limitingBox[coor], latticeVectors[coor]);
                }
            } catch (ArrayIndexOutOfBoundsException e) {
                System.err.println("Lx = " + Lx);
                System.err.println("Ly = " + Ly);
                System.err.println("n = " + n);
                System.err.println("k = " + k);
                System.err.println("j = " + j);
                System.err.println("i = " + i);
                System.err.println("l = " + l);
                System.err.println("sym_op = " + sym_op);
                e.printStackTrace();
            }
        }
        if (symVec[2]==1){
            loc = loc.subtract(latticeVectors[2].mapMultiply(2*k+1));
        }
        return loc;
    }

    public static int coordinateToInt(char r){
        switch(r) {
            case 'x': return 0;
            case 'y': return 1;
            case 'z': return 2;
        }
        // error
        return -1;
    }


    public double directSumNearestNeighbors(char[] interaction){
        RealVector origin = new ArrayRealVector(3); // constructs a vector of 3 zeros

        double sum = 0;
        int Lx=1; int Ly=1; int Lz=1;
        int[] neighbors=new int[]{1, 5, 19, 27};
        for (int i : neighbors) {
            RealVector ri = getLocation(i, new int[]{Lx, Ly});
//            System.out.println(ri.getNorm());
            sum += dipolar(coordinateToInt(interaction[2]), coordinateToInt(interaction[3]), origin, ri);
        }
        return sum;
    }


    public double directSum(int LxStart, int LxEnd, double heightTol, double radiusTol, char[] interaction){
        RealVector origin = new ArrayRealVector(3); // constructs a vector of 3 zeros
        double prevSum=0, LzSum=0;
        for (int Lx=LxStart; Lx<=LxEnd; Lx++){
            int Ly=Lx;
            double prevLzSum=0;
            LzSum=0;
            int i;
            for (i=1;; i++) {
                RealVector ri = getLocation(i, new int[]{Lx, Ly});
                if (ri.getSubVector(0, 2).getNorm() < Lx) {
                    if (interaction.length == 4) {
                        LzSum += dipolar(coordinateToInt(interaction[2]), coordinateToInt(interaction[3]), origin, ri);
                    } else if (interaction[4] == 'i' && interaction[5] == 'j'){
                        LzSum += dipolar(coordinateToInt(interaction[2]), coordinateToInt(interaction[3]), origin, ri)*dipolar(coordinateToInt(interaction[6]), coordinateToInt(interaction[7]), origin, ri);
                    } else if (interaction[4] == 'i' && interaction[5] == 'k'){
                        LzSum += calcPair(origin, ri, 0.05, 0.05, interaction, true);
                    } else if (interaction[4] == 'j' && interaction[5] == 'k'){
                        LzSum += calcPair(origin, ri, 0.05, 0.05, interaction, false);
                    }
//                    LzSum += dipolar(1, 2, origin, ri);//*dipolar(1, 0, origin, ri);
//                    LzSum += calcPair(origin, ri, 0.05, 0.005);
                }
                if (i % (num_in_cell * Lx * Ly * 8) == 0 && i / (num_in_cell * Lx * Ly * 8) > Lx) {    // a new "floor" is added at the top and bottom AND the height of the cylinder is larger than its basis radius
//                    System.out.println("Lx\t" + Lx + " Lz\t" + (i / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                    if ((i>100 && Math.abs(LzSum) < 1e-15) || Math.abs((LzSum - prevLzSum)/LzSum) < heightTol) {  // convergence
//                        System.out.println("Lz convergence! \t Lx\t" + Lx + " Lz\t" + (i / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                        break;
                    }
                    prevLzSum = LzSum;
                }
            }
            if ((Lx>3 && Math.abs(LzSum) < 1e-15) || Math.abs((LzSum - prevSum)/LzSum) < radiusTol){
//                System.out.println("Lx convergence!! \t Lx\t" + Lx + " Lz\t" + (i / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                break;
            }
            prevSum=LzSum;
        }
        if (Math.abs(LzSum) < 1E-14){
            return 0;
        }
        return LzSum;
    }

    public double calcPair(RealVector ri, RealVector rj, double heightTol, double radiusTol, char[] interaction, boolean bothWithRi){
        double prevSum=0, LzSum;
        for (int Lx = (int)(ri.subtract(rj).getSubVector(0,2).getNorm()+1);;Lx++){
            int Ly = Lx;
            double prevLzSum=0;
            LzSum=0;
            int k;
            for (k=1;;k++){
                RealVector rk = getLocation(k, new int[]{Lx, Ly});
                if (rk.getSubVector(0, 2).getNorm() < Lx && !ri.equals(rk) && !rj.equals(rk)){
//                    LzSum += dipolar(0, 2, ri,rk) * dipolar(0, 2, rk,rj) + dipolar(1, 2, ri,rk) * dipolar(1, 2, rk,rj);
//                    System.out.println(dipolar(1, 2, ri,rk) * dipolar(0, 2, rk,rj));
                    if (bothWithRi) {
                        LzSum += dipolar(coordinateToInt(interaction[2]), coordinateToInt(interaction[3]), ri, rk) * dipolar(coordinateToInt(interaction[6]),coordinateToInt(interaction[7]), ri, rj);
                    } else {
                        LzSum += dipolar(coordinateToInt(interaction[2]), coordinateToInt(interaction[3]), ri, rk) * dipolar(coordinateToInt(interaction[6]),coordinateToInt(interaction[7]), rk, rj);
                    }
                }
                if (k % (num_in_cell * Lx * Ly * 8) == 0 && k / (num_in_cell * Lx * Ly * 8) > Lx) {    // a new "floor" is added at the top and bottom AND the height of the cylinder is larger than its basis radius
                    if ((Math.abs(LzSum) < 1e-15) || Math.abs((LzSum - prevLzSum)/LzSum) < heightTol) {  // convergence
//                        System.out.println("Lz convergence! \t Lx\t" + Lx + " Lz\t" + (k / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                        break;
                    }
                    prevLzSum = LzSum;
                }
            }
            if ((Math.abs(LzSum) < 1e-15) || Math.abs((LzSum - prevSum)/LzSum) < radiusTol){
//                System.out.println("Lx convergence!! \t Lx\t" + Lx + " Lz\t" + (k / (num_in_cell * Lx * Ly * 8)) + " Lz_sum\t " + LzSum);
                break;
            }
            prevSum=LzSum;
        }
        return LzSum;
    }

    public double calcPairSphere(RealVector ri, RealVector rj, double radiusTol){
        // find center
        RealVector center = ri.add(rj).mapDivide(2.0);
        // find closest lattice point to center
        int[] cntrCoordinates = new int[3];
        for (int c=0;c< cntrCoordinates.length; c++) {
            cntrCoordinates[c] = (int) Math.floor(center.dotProduct(latticeVectors[c].unitVector()) / latticeVectors[c].getNorm());
        }
        center = new ArrayRealVector(3);
        for (int c=0;c< cntrCoordinates.length; c++) {
            center.combineToSelf(1, cntrCoordinates[c], latticeVectors[c]);
        }
        double prevSum=0, sum;
        for (int L = 2*(int)Math.ceil(ri.subtract(rj).getNorm());;L++){
            // L is the linear size of the bounding box
            int k;
            sum=0;
            for (k=0;k<num_in_cell*L*L*L*8;k++){
                RealVector rkRelToCntr = getLocation(k, new int[]{L, L});
                RealVector rk = center.add(rkRelToCntr);
                if (rkRelToCntr.getNorm() < L && !ri.equals(rk) && !rj.equals(rk)){
                    sum += dipolar(0, 2, ri,rk) * dipolar(0, 2, rk,rj) + dipolar(1, 2, ri,rk) * dipolar(1, 2, rk,rj);
                }
            }
            if (Math.abs((sum - prevSum)/sum) < radiusTol){ // convergence
//                System.out.println("Radius convergence!! \t L\t" + L + " sum\t " + sum);
                break;
            }
            prevSum=sum;
        }
        return sum;
    }

    public static void main(String[] args){
        //int Lx = Integer.parseInt(args[0]);
        CrystalSum cs = new CrystalSum();
//        final ArrayList<char[]> interactions = new ArrayList<char[]>() {{
//            add(new char[] {'i','j','x','x','i','j','y','y'});
//            add(new char[] {'i','j','x','y','i','j','y','y'});
//            add(new char[] {'i','j','y','y','i','j','z','z'});
//            add(new char[] {'i','j','y','y','i','k','y','y'});
//            add(new char[] {'i','j','z','z','i','k','y','y'});
//            add(new char[] {'i','j','y','y','i','k','z','z'});
//            add(new char[] {'i','j','z','z','i','k','z','z'});
//            add(new char[] {'i','j','x','y','j','k','x','y'});
//            add(new char[] {'i','j','y','y','j','k','x','y'});
//            add(new char[] {'i','j','x','x','j','k','y','y'});
//            add(new char[] {'i','j','x','y','j','k','y','y'});
//            add(new char[] {'i','j','y','y','j','k','y','y'});
//            add(new char[] {'i','j','z','z','j','k','y','y'});
//            add(new char[] {'i','j','x','z','j','k','y','z'});
//            add(new char[] {'i','j','y','z','j','k','y','z'});
//            add(new char[] {'i','j','y','y','j','k','z','z'});
//            add(new char[] {'i','j','z','z','j','k','z','z'});
//            add(new char[] {'i','j','x','x'});
//            add(new char[] {'i','j','x','y'});
//            add(new char[] {'i','j','y','x'});
//            add(new char[] {'i','j','y','y'});
//            add(new char[] {'i','j','y','z','i','j','y','z'});
//            add(new char[] {'i','j','x','x'});
//            add(new char[] {'i','j','x','y'});
//            add(new char[] {'i','j','x','z'});
//            add(new char[] {'i','j','y','x'});
//            add(new char[] {'i','j','y','y'});
//            add(new char[] {'i','j','y','z'});
//            add(new char[] {'i','j','z','x'});
//            add(new char[] {'i','j','z','y'});
//            add(new char[] {'i','j','z','z'});
//            add(new char[] {'i','j','z','x'});
//            add(new char[] {'i','j','z','z'});
//            add(new char[] {'i','j','y','z'});
//            add(new char[] {'i','j','z','y'});
//        }};
//        System.out.println(cs.directSum(Lx, Lx, 0.00000005, 0.0000005));
        final ArrayList<char[]> interactions = new ArrayList<char[]>();
        char[] coords = new char[]{'x','y','z'};
        for (int i=0;i<coords.length;i++){
            for (int j=0;j<coords.length;j++){
                for (int k=0;k<coords.length;k++){
                    for (int l=0;l<coords.length;l++){
                        interactions.add(new char[] {'i','j',coords[i],coords[j],'i','j',coords[k],coords[l]});
                    }
                }
            }
        }

        interactions.stream().forEach( interaction ->{
            if (interaction.length==8) {
//                System.out.println(String.format("V[%c,%c,%c,%c] V[%c,%c,%c,%c] -> factor (", interaction[0], interaction[1], interaction[2], interaction[3], interaction[4], interaction[5], interaction[6], interaction[7])
                System.out.print(String.format("B[%c,%c,%c,%c] = factor (", interaction[2], interaction[3], interaction[6], interaction[7])
                        + cs.directSum(3, 8, 0.00000000000005, 0.00000000000005, interaction) + "); ");
            } else {
//                System.out.println(String.format("V[%c,%c,%c,%c] -> Sqrt[factor] (", interaction[0], interaction[1], interaction[2], interaction[3])
                System.out.print(String.format("A[%c,%c] = Sqrt[factor] (", interaction[2], interaction[3])
//                System.out.println(String.format("V[%c,%c,%c,%c] \\[Delta]nn[i, j] -> Sqrt[factor] (", interaction[0], interaction[1], interaction[2], interaction[3])
                            + cs.directSum(3, 8, 0.000000000005, 0.000000000005, interaction) + "); ");
//                            + cs.directSumNearestNeighbors(interaction) + "), ");
            }
        });
//        System.out.println(cs.calcPair(new ArrayRealVector(new double[]{0,0,0}), new ArrayRealVector(new double[]{1,0,0}), 0.00000005, 0.0005));
//        System.out.println(cs.calcPairSphere(new ArrayRealVector(new double[]{0,0,0}), new ArrayRealVector(new double[]{0.5,0,0.25*cs.c}), 0.0000005));
    }
}
