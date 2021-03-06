// Version for students

package tudelft.cgv.volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import tudelft.cgv.gui.RaycastRendererPanel;
import tudelft.cgv.gui.TransferFunction2DEditor;
import tudelft.cgv.gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import tudelft.cgv.util.TFChangeListener;
import tudelft.cgv.util.VectorMath;
import tudelft.cgv.volume.GradientVolume;
import tudelft.cgv.volume.Volume;
import tudelft.cgv.volume.VoxelGradient;

import java.awt.Color;


/**
 *
 * @author michel
 *  Edit by AVilanova & Nicola Pezzotti
 * 
 * 
 * Main functions to implement the volume rendering
 */


//////////////////////////////////////////////////////////////////////
///////////////// CONTAINS FUNCTIONS TO BE IMPLEMENTED ///////////////
//////////////////////////////////////////////////////////////////////

public class RaycastRenderer extends Renderer implements TFChangeListener {


// attributes
    
    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunction2D tFunc2D;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;
    private String shadingType = "Default";
    private boolean toneShadingMode = false;
    private boolean isoMode = false;
    private float iso_value=95; 
    // This is a work around
    private float res_factor = 1.0f;
    private float max_res_factor=0.25f;
    private TFColor isoColor; 

    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE MODIFIED    /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    //Function that updates the "image" attribute (result of renderings)
    // using the slicing technique. 
    
    public void slicer(double[] viewMatrix) {
	
        // we start by clearing the image
        resetImage();

        // vector uVec and vVec define the view plane, 
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);

        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points. 
        for(int k=0;k<3;k++)
        {
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }

        // compute the volume center
        double[] volumeCenter = new double[3];
        computeVolumeCenter(volumeCenter);

        // Here will be stored the 3D coordinates of every pixel in the plane 
        double[] pixelCoord = new double[3];

        // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();

        int[] imageCenter = new int[2];
        // Center of the image/texture 
        imageCenter[0]= imageW/2;
        imageCenter[1]= imageH/2;
        
        // imageW/ image H contains the real width of the image we will use given the resolution. 
        //The resolution is generated once based on the maximum resolution.
        imageW = (int) (imageW*((max_res_factor/res_factor)));
        imageH = (int) (imageH*((max_res_factor/res_factor)));

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        // Color that will be used as a result 
        TFColor pixelColor = new TFColor();
        // Auxiliar color
        TFColor colorAux;

        // Contains the voxel value of interest
        int val;
        
        //Iterate on every pixel
        for (int j = imageCenter[1] - imageH/2; j < imageCenter[1] + imageH/2; j++) {
            for (int i =  imageCenter[0] - imageW/2; i <imageCenter[0] + imageW/2; i++) {
                
                // computes the pixelCoord which contains the 3D coordinates of the pixels (i,j)
                computePixelCoordinatesFloat(pixelCoord,volumeCenter,uVec,vVec,i,j);
                
                //we now have to get the value for the in the 3D volume for the pixel
                //we can use a nearest neighbor implementation like this:
                //val = volume.getVoxelNN(pixelCoord);
                
                //you have also the function getVoxelLinearInterpolated in Volume.java          
                val = (int) volume.getVoxelLinearInterpolate(pixelCoord);
                
                //you have to implement this function below to get the cubic interpolation
                //val = (int) volume.getVoxelTriCubicInterpolate(pixelCoord);
                
                
                // Map the intensity to a grey value by linear scaling
                 pixelColor.r = (val/max);
                 pixelColor.g = pixelColor.r;
                 pixelColor.b = pixelColor.r;

                // the following instruction makes intensity 0 completely transparent and the rest opaque
                // pixelColor.a = val > 0 ? 1.0 : 0.0;   
                
                // Alternatively, apply the transfer function to obtain a color using the tFunc attribute
                // colorAux= tFunc.getColor(val);
                // pixelColor.r=colorAux.r;pixelColor.g=colorAux.g;pixelColor.b=colorAux.b;pixelColor.a=colorAux.a; 
                // IMPORTANT: You can also simply use pixelColor = tFunc.getColor(val); However then you copy by reference and this means that if you change 
                // pixelColor you will be actually changing the transfer function So BE CAREFUL when you do this kind of assignments

                //BufferedImage/image/texture expects a pixel color packed as ARGB in an int
                //use the function computeImageColor to convert your double color in the range 0-1 to the format need by the image
                int pixelColor_i = computeImageColor(pixelColor.r,pixelColor.g,pixelColor.b,pixelColor.a);
                image.setRGB(i, j, pixelColor_i);
            }
        }
}
    

    
    //Do NOT modify this function
    //
    //Function that updates the "image" attribute using the MIP raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
    int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
    	//compute the increment and the number of samples
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        
        // Compute the number of times we need to sample
        double distance = VectorMath.distance(entryPoint, exitPoint);
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);

        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
       
        double maximum = 0;
        do {
            double value = volume.getVoxelLinearInterpolate(currentPos)/255.; 
            if (value > maximum) {
                maximum = value;
            }
            for (int i = 0; i < 3; i++) {
                currentPos[i] += increments[i];
            }
            nrSamples--;
        } while (nrSamples > 0);

        double alpha;
        double r, g, b;
        if (maximum > 0.0) { // if the maximum = 0 make the voxel transparent
            alpha = 1.0;
        } else {
            alpha = 0.0;
        }
        r = g = b = maximum;
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }
    
          
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    //Function that updates the "image" attribute using the Isosurface raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
   public int traceRayIso(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {

       double[] lightVector = new double[3];
       //We define the light vector as directed toward the view point (which is the source of the light)
       // another light vector would be possible
       VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);

       double[] increments = new double[3];
       VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep,
               rayVector[2] * sampleStep);
       // Implementation

       // Compute distance between entry and exit
       double totalDistance = VectorMath.distance(entryPoint, exitPoint);

       // Compute number of samples using distance and sample step length
       int totalSamples = (int)Math.floor(totalDistance/sampleStep) + 1;

       // Start from the entry point
       double [] currentPosition = new double [3];
       VectorMath.setVector(currentPosition, entryPoint[0], entryPoint[1], entryPoint[2]);

       //Initialization of the colors as floating point values
       double r, g, b;
       r = g = b = 0.0;
       double alpha = 0.0;
       TFColor voxelColor = new TFColor();
       // Now loop through all samples
       do{
           // Obtain the interpolation value at current position and compare it against iso_value
           // Using Linear Interpolation because tricubic interpolation is too slow
           // Tri cubic interpolation may be activated by un-commenting the command below
           double currentValue = volume.getVoxelLinearInterpolate(currentPosition);
           //double currentValue = volume.getVoxelTriCubicInterpolate(currentPosition);

           // Edit color values if currentValue is greater than iso_value
           if (currentValue > iso_value)
           {
               // Use bisection accuracy to fine tune current position
               currentPosition = bisection_accuracy(currentPosition, increments, sampleStep, currentValue, iso_value);
               r = isoColor.r; g = isoColor.g; b = isoColor.b; alpha = 1.0;
               break;
           }
        
        // Increment current position and decrement available samples
           for (int i = 0; i < 3; i++)
           {
               currentPosition[i] += increments[i];
           }
           totalSamples --;

       }while (totalSamples > 0);

       TFColor newColor;
       if((r > 0 || g > 0 || b > 0) && alpha > 0) {
           if (shadingMode) {
               voxelColor = new TFColor(r,g,b,alpha);
               newColor = computePhongShading(voxelColor, gradients.getGradient(currentPosition), lightVector, rayVector, "Default");
               r=newColor.r;
               g=newColor.g;
               b=newColor.b;
               alpha = newColor.a;

           }
       }

       // isoColor contains the isosurface color from the interface
       //r = isoColor.r;g = isoColor.g;b =isoColor.b;alpha =1.0;
       //computes the color
       int color = computeImageColor(r,g,b,alpha);
       return color;
    }
   
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    
    // Given the current sample position, increment vector of the sample (vector from previous sample to current sample) and sample Step. 
   // Previous sample value and current sample value, isovalue value
    // The function should search for a position where the iso_value passes that it is more precise.
   public double[]  bisection_accuracy (double[] currentPos, double[] increments,double sampleStep, double value, double iso_value)
   {
       //Middle value for each iteration
       float midValue;

       //Flag to check if the value has been found
       boolean found = false;

       //Convergence value as a check to end iterating
       double convergence = 0;

       //previous position
       double[] previousPos = {currentPos[0] - increments[0], currentPos[1] - increments[1], currentPos[2] - increments[2]};

       do {
           //calculate the mid value and position between previous and current position
           double[] midPos = {(currentPos[0] + previousPos[0]) / 2, (currentPos[1] + previousPos[1]) / 2, (currentPos[2] + previousPos[2]) / 2};
           midValue = volume.getVoxelLinearInterpolate(midPos);

           //Calculate the difference in position magnitude as a measure for convergence
           convergence=Math.sqrt(Math.pow(currentPos[0], 2) + Math.pow(currentPos[1], 2) + Math.pow(currentPos[2], 2))-Math.sqrt(Math.pow(previousPos[0], 2) + Math.pow(previousPos[1], 2) + Math.pow(previousPos[2], 2));
           //Check if midPos is greater than, less than or  equal to iso value and update accordingly
           if(midValue == iso_value) {
               currentPos = midPos;
               found = true;
           }
           else if (midValue > iso_value) {
               currentPos = midPos;
           }
           else {
               previousPos = midPos;
           }

       } while (found && convergence > 0.001); 
       
       return currentPos;
   }
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    //Function that updates the "image" attribute using the compositing// accumulatted raycasting
    //It returns the color assigned to a ray/pixel given it's starting point (entryPoint) and the direction of the ray(rayVector).
    // exitPoint is the last point.
    //ray must be sampled with a distance defined by the sampleStep
   
    public int traceRayComposite(double[] entryPoint, double[] exitPoint, double[] rayVector, double sampleStep) {
        double[] lightVector = new double[3];
        
        //the light vector is directed toward the view point (which is the source of the light)
        // another light vector would be possible 
        VectorMath.setVector(lightVector, rayVector[0], rayVector[1], rayVector[2]);
        
        // Compute the number of times we need to sample
        int nrSamples = 1 + (int) Math.floor(VectorMath.distance(entryPoint, exitPoint) / sampleStep);
        
        // Define increments
        double[] increments = new double[3];
        VectorMath.setVector(increments, rayVector[0] * sampleStep, rayVector[1] * sampleStep, rayVector[2] * sampleStep);
        
        //the current position is initialized as the entry point
        double[] currentPos = new double[3];
        VectorMath.setVector(currentPos, entryPoint[0], entryPoint[1], entryPoint[2]);
        
        //Initialization of the colors as floating point values
        double r, g, b;
        r = g = b = 0.0;
        double alpha = 0.0; // alpha is opacity here
        double opacity = 0;
        
        TFColor voxelColor = new TFColor();
        TFColor compositeColor = new TFColor();
 
        if (compositingMode) {
            // 1D transfer function 
            compositeColor = getCompositeColor(currentPos, lightVector, nrSamples, rayVector);
        }
                    
        if (tf2dMode) {
            compositeColor = getCompositeColor2D(currentPos, lightVector, nrSamples, rayVector);
        }
        
        // copy the results.
        voxelColor.r = compositeColor.r;
        voxelColor.g = compositeColor.g;
        voxelColor.b = compositeColor.b;
        // only give opacity if it has some color.
        if(voxelColor.r > 0 || voxelColor.g > 0 || voxelColor.b > 0) {
            opacity = compositeColor.a;
        }
        
        r = voxelColor.r ;
        g = voxelColor.g ;
        b = voxelColor.b;
        alpha = opacity ;
            
        //computes the color
        int color = computeImageColor(r,g,b,alpha);
        return color;
    }

    
    // Get the composite color from raycasting front to back.
    // Get the composite color from raycasting front to back.
    public TFColor getCompositeColor(double[] currentPos, double[] lightVector, int nrSamples, double[] rayVector){

        int value = (int) volume.getVoxelLinearInterpolate(currentPos);
        
        // Use color from tFunc widget.
        TFColor currColor = this.tFunc.getColor(value);
        double currentOpacity = currColor.a;
        
        // set an opacity threshold when computation can stop.
        double opacity_threshold = 1.5;
        
        // End the ray when opacity is too high.
        // Get just this sample's color value * opacity.
        if (nrSamples < 0 || currentOpacity > opacity_threshold) {
            currColor.r *= currentOpacity;
            currColor.g *= currentOpacity;
            currColor.b *= currentOpacity;
            return currColor;
        }
        
        if(toneShadingMode) {
            // Re-compute currColor to avoid pass by reference from other modules
            currColor = this.tFunc.getColor(value);
            currColor = computeToneShading(currColor, this.gradients.getGradient(currentPos), lightVector, rayVector);
        } else if (shadingMode) {
            // Re-compute currColor to avoid pass by reference from other modules
            currColor = this.tFunc.getColor(value);
            currColor = computePhongShading(currColor, this.gradients.getGradient(currentPos), 
                    lightVector, rayVector, shadingType);
            currentOpacity = currColor.a;
        }
        
        
        // Go to the next position.
        for (int i = 0; i < 3; i++) {
            currentPos[i] += lightVector[i];
        }
        nrSamples--;

        // Get the accumulated color of the sample before.
        TFColor previousColor = getCompositeColor(currentPos, lightVector, nrSamples, rayVector);
        
        // Update the color with behind sample's color.
        // The amount of previous color that can shine through, is  this color's transparancy.
        TFColor newColor = new TFColor();
        newColor.r = currentOpacity * currColor.r + (1 - currentOpacity) * previousColor.r;
        newColor.g = currentOpacity * currColor.g + (1 - currentOpacity) * previousColor.g;
        newColor.b = currentOpacity * currColor.b + (1 - currentOpacity) * previousColor.b;
        newColor.a = currentOpacity + (1 - currentOpacity) * previousColor.a;
        
        return newColor;
    }
            
    public TFColor getCompositeColor2D(double[] currentPos, double[] lightVector, int nrSamples, double[] rayVector){
        
        if (nrSamples < 0) {
            return new TFColor(0,0,0,0);
        }
        
        // Voxel value, gradient, and gradient magnitude for 2D transfer function
        int value = (int) volume.getVoxelLinearInterpolate(currentPos);
        VoxelGradient currentGradient = gradients.getGradient(currentPos);
        double mag = currentGradient.mag;

        // Go to the next position.
        for (int i = 0; i < 3; i++) {
            currentPos[i] += lightVector[i];
        }
        nrSamples--;
        
        // Use color from tFund2D widget and opacity from computeOpacity2D
        TFColor currColor = this.tFunc2D.color;
        double currentOpacity = this.computeOpacity2DTF(tFunc2D.baseIntensity, tFunc2D.radius, 
                value, "Default", rayVector, currentGradient)*currColor.a;
        // Get the accumulated color of the sample before.
        TFColor previousColor = getCompositeColor2D(currentPos, lightVector, nrSamples, rayVector);
           
        // Add shading if necessary
        if (shadingMode) {
            currColor = computePhongShading(currColor, gradients.getGradient(currentPos), lightVector, rayVector, "Default");
            // Implementing slightly tuned versions of improved opacity computation to render finer details
            if (shadingType == "Default"){
                currentOpacity = Math.min(this.computeOpacity2DTF(tFunc2D.baseIntensity, tFunc2D.radius, value, 
                    shadingType, rayVector, currentGradient)*currColor.a, 1.5);
            }
            else if (shadingType == "Surface"){
                currentOpacity = Math.min(this.computeOpacity2DTF(tFunc2D.baseIntensity, tFunc2D.radius, value, 
                    shadingType, rayVector, currentGradient)*currColor.a, 0.3);
            }
            else if (shadingType == "Boundary"){
                currentOpacity = Math.min(this.computeOpacity2DTF(tFunc2D.baseIntensity, tFunc2D.radius, value, 
                    shadingType, rayVector, currentGradient)*currColor.a, 0.3);
            }
        }
        
        // Update the color with behind sample's color.
        // The amount of previous color that can shine through, is  this color's transparancy.
        TFColor newColor = new TFColor(0,0,0,0);
        newColor.r = currentOpacity * currColor.r + (1 - currentOpacity) * previousColor.r;
        newColor.g = currentOpacity * currColor.g + (1 - currentOpacity) * previousColor.g;
        newColor.b = currentOpacity * currColor.b + (1 - currentOpacity) * previousColor.b;
        // Compute accumulated opacity using formula
        newColor.a = currentOpacity + (1 - currentOpacity) * previousColor.a;
        
        return newColor;
    }
    
    // Method to apply tone shading.
    public TFColor computeToneShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector,
            double[] rayVector) {
        // define cool and warm colors.
        TFColor cool = new TFColor(0.4, 0.0, 1.0, 1.0);
        TFColor warm = new TFColor(1.0, 1.0, 0.0, 1.0);
        
        // weight of TF color compared to tone.
        double weight = 0.1;
        
        // Normalise Gradient Vector
        double gradientNorm[] = new double[3];
        if (gradient.mag != 0){
            VectorMath.setVector(gradientNorm, gradient.x/gradient.mag, gradient.y/gradient.mag, gradient.z/gradient.mag);
        } else{
            return voxel_color;
        }
        
        // Normalise Gradient Vector
        double lightNorm[] = new double[3];
        if (VectorMath.length(lightVector) != 0){
            VectorMath.setVector(lightNorm, lightVector[0]/VectorMath.length(lightVector), lightVector[1]/VectorMath.length(lightVector), lightVector[2]/VectorMath.length(lightVector));
        } else {
            return voxel_color;
        }

        // Calculate cos theta using gradient (normal) vector
        double cosTheta = VectorMath.dotproduct(lightNorm, gradientNorm);
        
        double coolColor[] = new double[3];
        double warmColor[] = new double[3];
        
        // Compute the combined color.
        VectorMath.setVector(coolColor, cool.r + weight * voxel_color.r, 
                cool.g + weight * voxel_color.g, cool.b + weight * voxel_color.b);
        VectorMath.setVector(warmColor, warm.r + weight * voxel_color.r, 
                warm.g + weight * voxel_color.g, warm.b + weight * voxel_color.b);
        
        // Define the cool to warm ratio.
        double coolComp = (1 + cosTheta) / 2;
        double warmComp = 1 - ((1 + cosTheta) / 2);
        
        // Apply the ratio to the warm and cool components.
        double r = coolComp * coolColor[0] + warmComp * warmColor[0];
        double g = coolComp * coolColor[1] + warmComp * warmColor[1];
        double b = coolComp * coolColor[2] + warmComp * warmColor[2];
        
        TFColor color = new TFColor(r,g,b,voxel_color.a);
        return color;
    }
    
    //////////////////////////////////////////////////////////////////////
    ///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Compute Phong Shading given the voxel color (material color), the gradient, the light vector and view vector 
    public TFColor computePhongShading(TFColor voxel_color, VoxelGradient gradient, double[] lightVector, 
            double[] rayVector, String shadingType) {

        // Implementation
        // Define ambient, diffuse, and specular constants
        double ka = 0.1;
        double kd = 0.7;
        double ks = 0.2;
        double alpha = 100.0;

        double gradientNorm[] = new double[3];
        // Normalise Gradient Vector
        if (gradient.mag != 0){
            VectorMath.setVector(gradientNorm, gradient.x/gradient.mag, gradient.y/gradient.mag, gradient.z/gradient.mag);
        }
        // Return voxel color is no gradient at current position
        else{
            return voxel_color;
        }

        // Calculate cos theta using gradient (normal) vector
        double cosTheta = VectorMath.dotproduct(lightVector, gradientNorm);

        // Calculate cos phi
        // Scale up light vector by twice the scalar alpha
        double n_dot_L = VectorMath.dotproduct(gradientNorm, lightVector);

        double [] scaledNormalVec = {2*n_dot_L*gradientNorm[0], 2*n_dot_L*gradientNorm[1], 2*n_dot_L*gradientNorm[2]};

        double [] vectorSubstraction = {scaledNormalVec[0]-lightVector[0], scaledNormalVec[1]-lightVector[1],
                scaledNormalVec[2]-lightVector[2]};
        // Dot product of resultant and ray vector gives cos phi
        double cosPhi = VectorMath.dotproduct(rayVector, vectorSubstraction);
        
        // Compute ambient, diffused, and specular color components
        double ciAmbient[] = new double[3];
        VectorMath.setVector(ciAmbient, voxel_color.r,voxel_color.g,voxel_color.b);

        double ciDiffused[] = new double[3];
        VectorMath.setVector(ciDiffused, voxel_color.r,voxel_color.g,voxel_color.b);

        double ciSpecular[] = new double[3];
        VectorMath.setVector(ciSpecular, 1.0,1.0,1.0);

        double ambientColor[] = new double[3];
        double diffuseColor[] = new double[3];
        double specularColor[] = new double[3];

        for(int i= 0;i<3;i++){ambientColor[i] = Math.max(0,ka*ciAmbient[i]);}
        for(int i= 0;i<3;i++){diffuseColor[i] = Math.max(0,kd*cosTheta*ciDiffused[i]);}
        for(int i= 0;i<3;i++){specularColor[i] = Math.max(0,ks*Math.pow(cosPhi, alpha)*ciSpecular[i]);}
        
        // Compose final color by combining all color components
        double r = ambientColor[0] + diffuseColor[0] + specularColor[0];
        if (r < 0.0) {r = 0.0;}
        if (r > 1.0) {r = 1.0;}
        double g = ambientColor[1] + diffuseColor[1] + specularColor[1];
        if (g < 0.0) {g = 0.0;}
        if (g > 1.0) {g = 1.0;}
        double b = ambientColor[2] + diffuseColor[2] + specularColor[2];
        if (b < 0.0) {b = 0.0;}
        if (b > 1.0) {b = 1.0;}
        
        // Compute improved opacity based on the user input
        double opacity = improveOpacity(voxel_color.a, gradient, rayVector, shadingType);
        TFColor color = new TFColor(r,g,b,opacity);
        return color;
    }

    public double improveOpacity(double currentOpacity, VoxelGradient gradient, double[] rayVector, String shadingType)
    {
        double gradientNorm[] = new double[3];
        if (gradient.mag != 0){
            VectorMath.setVector(gradientNorm, gradient.x/gradient.mag, gradient.y/gradient.mag, gradient.z/gradient.mag);
        } else{
            return currentOpacity;
        }
        // Perform silhouette rendering or not depending on the input
        double improvedOpacity = 0.0;
        if (shadingType == "Surface"){
            // Define constants for silhouette opacity
            double ksc = 0.0;
            double kss = 50;
            double kse = 0.6;
            double n_dot_R = VectorMath.dotproduct(gradientNorm, rayVector);
            // Calculate improved opacity
            improvedOpacity = Math.min(currentOpacity*(ksc + kss*Math.pow((1 - Math.abs(n_dot_R)), kse)), 1.2);
        } else if (shadingType == "Boundary"){
            // Define constants for gradient opacity
            double kgc = 0.0;
            double kgs = 0.05;
            double kge = 1.2;
            // Parameters for gradient opacity tuned
            improvedOpacity = Math.min(currentOpacity*(kgc + kgs*Math.pow(gradient.mag, kge)), 1.5);
        } else if (shadingType == "Default"){
          // Return current opacity value if default option selected
          improvedOpacity = currentOpacity;
        }
      
        return improvedOpacity;
    }

    //////////////////////////////////////////////////////////////////////
    ///////////////// LIMITED MODIFICATION IS NEEDED /////////////////////
    ////////////////////////////////////////////////////////////////////// 
    // Implements the basic tracing of rays trough the image and given the
    // camera transformation
    // It calls the functions depending on the raycasting mode
  
    public void raycast(double[] viewMatrix) {

    	//data allocation
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        
        // increment in the pixel domain in pixel units
        int increment = 1;
        // sample step in voxel units
        int sampleStep = 1;
        // reset the image to black
        resetImage();
        
        // vector uVec and vVec define the view plane, 
        // perpendicular to the view vector viewVec which is going from the view point towards the object
        // uVec contains the up vector of the camera in world coordinates (image vertical)
        // vVec contains the horizontal vector in world coordinates (image horizontal)
        getViewPlaneVectors(viewMatrix,viewVec,uVec,vVec);
        
        
        // The result of the visualization is saved in an image(texture)
        // we update the vector according to the resolution factor
        // If the resolution is 0.25 we will sample 4 times more points. 
        for(int k=0;k<3;k++)
        {
            uVec[k]=res_factor*uVec[k];
            vVec[k]=res_factor*vVec[k];
        }
        
       // We get the size of the image/texture we will be puting the result of the 
        // volume rendering operation.
        int imageW=image.getWidth();
        int imageH=image.getHeight();

        int[] imageCenter = new int[2];
        // Center of the image/texture 
        imageCenter[0]= imageW/2;
        imageCenter[1]= imageH/2;
        
        // imageW/ image H contains the real width of the image we will use given the resolution. 
        //The resolution is generated once based on the maximum resolution.
        imageW = (int) (imageW*((max_res_factor/res_factor)));
        imageH = (int) (imageH*((max_res_factor/res_factor)));
        
        //The rayVector is pointing towards the scene
        double[] rayVector = new double[3];
        rayVector[0]=-viewVec[0];rayVector[1]=-viewVec[1];rayVector[2]=-viewVec[2];
             
        // compute the volume center
        double[] volumeCenter = new double[3];
        computeVolumeCenter(volumeCenter);

        
        // ray computation for each pixel
        for (int j = imageCenter[1] - imageH/2; j < imageCenter[1] + imageH/2; j += increment) {
            for (int i =  imageCenter[0] - imageW/2; i <imageCenter[0] + imageW/2; i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
            	computePixelCoordinatesBehindFloat(pixelCoord,viewVec,uVec,vVec,i,j);
            	// compute the entry and exit point of the ray
                computeEntryAndExit(pixelCoord, rayVector, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int val = 0;
                    if (compositingMode || tf2dMode) {
                        val = traceRayComposite(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (mipMode) {
                        val = traceRayMIP(entryPoint, exitPoint, rayVector, sampleStep);
                    } else if (isoMode){
                        val= traceRayIso(entryPoint,exitPoint,rayVector, sampleStep);
                    }
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, val);
                        }
                    }
                }

            }
        }
    }


//////////////////////////////////////////////////////////////////////
///////////////// FUNCTION TO BE IMPLEMENTED /////////////////////////
////////////////////////////////////////////////////////////////////// 
// Compute the opacity based on the value of the pixel and the values of the
// triangle widget tFunc2D contains the values of the baseintensity and radius
// tFunc2D.baseIntensity, tFunc2D.radius they are in image intensity units

public double computeOpacity2DTF(double material_value, double material_r,
        double voxelValue, String shadingType, double [] rayVector, VoxelGradient gradient) {

    //init opacity with 0
    double opacity = 0.0;

    //maximum magnitude
    double maxmag = this.gradients.getMaxGradientMagnitude();
    //angle of the widget
    double angle = Math.atan(material_r/maxmag);
    
    // Compute gradient magnitude
    double gradMagnitude = gradient.mag;
    //angle of current voxel with respect to base intensity center
    double voxelRad = Math.abs(voxelValue-material_value);
    double voxelAngle = Math.atan(voxelRad/gradMagnitude);
    
    double voxelBaseRadius = (gradMagnitude/maxmag)*material_r*0.5;

    //if the voxel is inside the widget, give it an opacity
    if(voxelAngle < angle){
        //the factor between the voxel radius and the total distance of the horizontal line to the diagonal along the voxel is used as a ramp
        opacity = tFunc2D.color.a;
        //opacity = (voxelRad/voxelBaseRadius)*tFunc2D.color.a;
    }   
    double improvedOpacity = improveOpacity(opacity, gradient, rayVector, shadingType);
    return improvedOpacity;
}  

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
    
    //Do NOT modify this function
    int computeImageColor(double r, double g, double b, double a){
		int c_alpha = 	a <= 1.0 ? (int) Math.floor(a * 255) : 255;
        int c_red = 	r <= 1.0 ? (int) Math.floor(r * 255) : 255;
        int c_green = 	g <= 1.0 ? (int) Math.floor(g * 255) : 255;
        int c_blue = 	b <= 1.0 ? (int) Math.floor(b * 255) : 255;
        int pixelColor = getColorInteger(c_red,c_green,c_blue,c_alpha);
        return pixelColor;
	}
    //Do NOT modify this function    
    public void resetImage(){
    	for (int j = 0; j < image.getHeight(); j++) {
	        for (int i = 0; i < image.getWidth(); i++) {
	            image.setRGB(i, j, 0);
	        }
	    }
    }
   
    //Do NOT modify this function
    void computeIncrementsB2F(double[] increments, double[] rayVector, double sampleStep) {
        // we compute a back to front compositing so we start increments in the oposite direction than the pixel ray
    	VectorMath.setVector(increments, -rayVector[0] * sampleStep, -rayVector[1] * sampleStep, -rayVector[2] * sampleStep);
    }
    
    //used by the slicer
    //Do NOT modify this function
    void getViewPlaneVectors(double[] viewMatrix, double viewVec[], double uVec[], double vVec[]) {
            VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
	    VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
	    VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
	}
    
    //used by the slicer	
    //Do NOT modify this function
    void computeVolumeCenter(double volumeCenter[]) {
	VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
    }
	
    //used by the slicer
    //Do NOT modify this function
    void computePixelCoordinatesFloat(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], float i, float j) {
        // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
            float imageCenter = image.getWidth()/2;
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
    }
        
    //Do NOT modify this function
    void computePixelCoordinates(double pixelCoord[], double volumeCenter[], double uVec[], double vVec[], int i, int j) {
        // Coordinates of a plane centered at the center of the volume (volumeCenter and oriented according to the plane defined by uVec and vVec
            int imageCenter = image.getWidth()/2;
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + volumeCenter[0];
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + volumeCenter[1];
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + volumeCenter[2];
    }
    //Do NOT modify this function
    void computePixelCoordinatesBehindFloat(double pixelCoord[], double viewVec[], double uVec[], double vVec[], float i, float j) {
            int imageCenter = image.getWidth()/2;
            // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
            // to the diaganal to make sure I am far away enough.

            double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
    }
    //Do NOT modify this function
    void computePixelCoordinatesBehind(double pixelCoord[], double viewVec[], double uVec[], double vVec[], int i, int j) {
            int imageCenter = image.getWidth()/2;
            // Pixel coordinate is calculate having the center (0,0) of the view plane aligned with the center of the volume and moved a distance equivalent
            // to the diaganal to make sure I am far away enough.

            double diagonal = Math.sqrt((volume.getDimX()*volume.getDimX())+(volume.getDimY()*volume.getDimY())+ (volume.getDimZ()*volume.getDimZ()))/2;               
            pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * diagonal + volume.getDimX() / 2.0;
            pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * diagonal + volume.getDimY() / 2.0;
            pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * diagonal + volume.getDimZ() / 2.0;
    }
	
    //Do NOT modify this function
    public int getColorInteger(int c_red, int c_green, int c_blue, int c_alpha) {
    	int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    	return pixelColor;
    } 
    //Do NOT modify this function
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
        isoColor = new TFColor();
        isoColor.r=1.0;isoColor.g=1.0;isoColor.b=0.0;isoColor.a=1.0;
    }
    
     
    //Do NOT modify this function
    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ())* (1/max_res_factor));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
           
        // Initialize transferfunction 
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tFunc2D= new TransferFunction2D((short) (volume.getMaximum() / 2), 0.2*volume.getMaximum());
        tfEditor2D = new TransferFunction2DEditor(tFunc2D,volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }
    //Do NOT modify this function
    public RaycastRendererPanel getPanel() {
        return panel;
    }

    //Do NOT modify this function
    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    //Do NOT modify this function
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    //Do NOT modify this function
    public void setShadingMode(boolean mode, String inputShadingType) {
        setTotalShadingMode(mode, !mode, inputShadingType);
    }
    
    public void setToneShadingMode(boolean mode, String inputShadingType) {
        setTotalShadingMode(!mode, mode, inputShadingType);
    }
    
    public void setTotalShadingMode(boolean volumeMode, boolean toneMode, String inputShadingType){
        shadingMode = volumeMode;
        shadingType = inputShadingType;
        toneShadingMode = toneMode;
        changed();
    }

    //Do NOT modify this function
    public void setMIPMode() {
        setMode(false, true, false, false,false);
    }
    //Do NOT modify this function
    public void setSlicerMode() {
        setMode(true, false, false, false,false);
    }
    //Do NOT modify this function
    public void setCompositingMode() {
        setMode(false, false, true, false,false);
    }
    //Do NOT modify this function
    public void setTF2DMode() {
        setMode(false, false, false, true, false);
    }
    //Do NOT modify this function
    public void setIsoSurfaceMode(){
        setMode(false, false, false, false, true);
     }
    //Do NOT modify this function
    public void setIsoValue(float pIsoValue){
         iso_value = pIsoValue;
         if (isoMode){
             changed();
         }
             
     }
    //Do NOT modify this function
    public void setResFactor(int value) {
         float newRes= 1.0f/value;
         if (res_factor != newRes)
         {
             res_factor=newRes;
             if (volume != null) changed();
         }
     }
     
   //Do NOT modify this function
   public void setIsoColor(TFColor newColor)
     {
         this.isoColor.r=newColor.r;
         this.isoColor.g=newColor.g;
         this.isoColor.b=newColor.b;
         if ((volume!=null) && (this.isoMode)) changed();
     }
     
    //Do NOT modify this function
     public float getIsoValue(){
         return iso_value;
     }
    //Do NOT modify this function
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d, boolean iso) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        isoMode = iso;
        changed();
    }
    //Do NOT modify this function
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }
    //Do NOT modify this function
    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }
    //Do NOT modify this function
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) < 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
    
     
    
    //Do NOT modify this function
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }

    //Do NOT modify this function
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();


    }
    //Do NOT modify this function
    @Override
    public void visualize(GL2 gl) {

        double[] viewMatrix = new double[4 * 4];
        
        if (volume == null) {
            return;
        }
        	
         drawBoundingBox(gl);

        
     //    gl.glGetDoublev(GL2.GL_PROJECTION_MATRIX,viewMatrix,0);

         gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);
                      
         
        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);
        
        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        //gl.glEnable(GL.GL_BLEND);
        //gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        
        double halfWidth = res_factor*image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST);
        gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(texture.getImageTexCoords().left(), texture.getImageTexCoords().top());
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(texture.getImageTexCoords().left(), texture.getImageTexCoords().bottom());
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(texture.getImageTexCoords().right(), texture.getImageTexCoords().bottom());
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(texture.getImageTexCoords().right(), texture.getImageTexCoords().top());
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }


    }
    
    //Do NOT modify this function
    public BufferedImage image;

    //Do NOT modify this function
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }

}
