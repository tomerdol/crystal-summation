import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

public class CrystalSum {


    public static double dipolar(int comp1, int comp2, RealVector ri, RealVector rj){
        double norm = ri.subtract(rj).getNorm();
        if (comp1==comp2){
            return (norm*norm - 3*Math.pow(ri.subtract(rj).getEntry(comp1), 2))/Math.pow(norm, 5);
        } else {
            return (-3*ri.subtract(rj).getEntry(comp1)*ri.subtract(rj).getEntry(comp2))/Math.pow(norm, 5);
        }
    }

    public static RealVector getLocation(int n, int Lx, int Ly, int Lz, int num_in_cell, RealVector[] latticeVectors, double[][] basis){
        int k = ((n/num_in_cell)/Lx)/Ly/8%Lz;
        int sym_op = ((n/num_in_cell)/Lx)/Ly%8;
        int j = ((n/num_in_cell)/Lx)%Ly;
        int i = (n/num_in_cell)%Lx;
        int l = n%num_in_cell;

        int[] symVec = new int[]{sym_op%2,sym_op/2%2,sym_op/4%2};

        int[] locCoordinates = new int[]{i, j, k};
        RealVector loc = new ArrayRealVector(3);
        for (int coor=0;coor<latticeVectors.length;coor++){
            loc=loc.add(locCoordinates[coor], latticeVectors[coor]);	// add primitive lattice vectors according to (i,j,k)
            loc=loc.add(basis[l][coor],latticeVectors[coor]);	// add fraction of primitive lattice vector according to basis index l
        }
        RealVector position = latticeVectors[0]*i +
    }
    sym_vec=np.array(())
            pos = (lattice_vectors[0]*i+lattice_vectors[1]*j+lattice_vectors[2]*k+basis_vectors[l]
            - sym_vec[0]*lattice_vectors[0]*Lx - sym_vec[1]*lattice_vectors[1]*Ly)
    if sym_vec[2]==1:
    pos = pos - (2*k+1)*lattice_vectors[2]
            return pos

    public static void main(String[] args){
        double a=1;
        double c=2.077294686;
        int num_in_cell=4;
        RealVector[] latticeVectors = new RealVector[3];
        latticeVectors[0] = new ArrayRealVector(new double[]{a, 0, 0});
        latticeVectors[1] = new ArrayRealVector(new double[]{0, a, 0});
        latticeVectors[2] = new ArrayRealVector(new double[]{0, 0, c});
        double[][] basis = {
            {0,0,0},
            {0.5,0,0.25},
            {0.5,0.5,0.5},
            {0,0.5,0.75}
        };

        RealVector origin = new ArrayRealVector(3); // constructs a vector of 3 zeros

    }
}
