/*
 * Anna Vilanova: This class implements the basic interaction with the 2D Transfer
 * Function
 * 
 * YOU MIGHT WANT TO MODIFY IT FOR THE EXTENSIONS
 *
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package tudelft.cgv.volvis;

/**
 *
 * @author Anna Vilanova
 */
public class TransferFunction2D {
    // base of the intensity of the triangle
    public short baseIntensity;
    // width of half of the triangle in intensity value units 
    public double radius;
    public TFColor color;
        

        public TransferFunction2D(short base, double r) {
            this.baseIntensity = base;
            this.radius = r;
            this.color = new TFColor(0.0, 204.0/255.0, 153.0/255.0, 0.3);
        }
        
        public void SetBaseRadius(short base, double r)
        {   
            this.baseIntensity = base;
            this.radius = r;
        }
}
