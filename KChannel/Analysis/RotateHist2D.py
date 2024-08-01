import ROOT
import numpy as np
from array import array

# Function to rotate a 2D histogram
def RotateHistogram(hist, angle_degrees):
    # Clone the original histogram
    rotated_hist = hist.Clone(hist.GetName() + "_rotated")
    
    # Convert angle from degrees to radians
    angle_radians = ROOT.TMath.DegToRad() * angle_degrees
    
    # Rotation matrix
    rotation_matrix = ROOT.TMatrixD(2, 2)
    rotation_matrix[0][0] = ROOT.TMath.Cos(angle_radians)
    rotation_matrix[0][1] = -ROOT.TMath.Sin(angle_radians)
    rotation_matrix[1][0] = ROOT.TMath.Sin(angle_radians)
    rotation_matrix[1][1] = ROOT.TMath.Cos(angle_radians)
    
    # Get the bin contents and errors of the original histogram
    for binx in range(1, hist.GetNbinsX() + 1):
        for biny in range(1, hist.GetNbinsY() + 1):
            content = hist.GetBinContent(binx, biny)
            error = hist.GetBinError(binx, biny)
            
            # Get the center coordinates of the bin
            x_center = hist.GetXaxis().GetBinCenter(binx)
            y_center = hist.GetYaxis().GetBinCenter(biny)
            
            # Apply rotation to the coordinates
            x_rotated = rotation_matrix[0][0] * x_center + rotation_matrix[0][1] * y_center
            y_rotated = rotation_matrix[1][0] * x_center + rotation_matrix[1][1] * y_center
            
            # Set the bin content and error in the rotated histogram
            rotated_hist.Fill(x_rotated, y_rotated, content)
            rotated_hist.SetBinError(rotated_hist.FindBin(x_rotated, y_rotated), error)
    
    return rotated_hist

# Function to calculate ellipse parameters
def calculate_ellipse_parameters(hist, angle_degrees, confidence_level=0.99):
    # Convert angle from degrees to radians
    angle_radians = ROOT.TMath.DegToRad() * angle_degrees

    # Rotate the covariance matrix
    covariance_matrix = np.array([[hist.GetCovariance(1, 1), hist.GetCovariance(1, 2)],
                                  [hist.GetCovariance(2, 1), hist.GetCovariance(2, 2)]])
    
    rotation_matrix = np.array([[ROOT.TMath.Cos(angle_radians), -ROOT.TMath.Sin(angle_radians)],
                                [ROOT.TMath.Sin(angle_radians), ROOT.TMath.Cos(angle_radians)]])
    
    rotated_covariance_matrix = np.dot(rotation_matrix, np.dot(covariance_matrix, rotation_matrix.T))

    # Perform Singular Value Decomposition (SVD) to get principal components
    singular_values, principal_components = np.linalg.eig(rotated_covariance_matrix)
    ellipse_width = np.sqrt(confidence_level * singular_values[0])
    ellipse_height = np.sqrt(confidence_level * singular_values[1])

    # Calculate the center of the ellipse
    center_x = hist.GetMean(1)
    center_y = hist.GetMean(2)

    # Draw the ellipse
    ellipse = ROOT.TEllipse(center_x, center_y, ellipse_width, ellipse_height, 0, 360, angle_degrees)
    ellipse.SetLineColor(ROOT.kRed)
    ellipse.SetLineWidth(2)
    ellipse.SetFillStyle(0)
    ellipse.Draw()

# Example usage
def RotatedHist():
    # Create a unique name for the histogram
    hist_name = "hist_" + str(ROOT.TUUID().AsString())

    # Create a 2D histogram
    hist = ROOT.TH2D(hist_name, "2D Histogram", 50, -5, 5, 50, -5, 5)
    
    # Fill histogram with random Gaussian-distributed data
    rng = ROOT.TRandom3()
    for i in range(10000):
        x = rng.Gaus(0, 1)
        y = rng.Gaus(0, 1)
        hist.Fill(x, y)
    
    # Rotate the histogram
    rotated_hist = RotateHistogram(hist, 45) # Rotate by 45 degrees
    
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", "Canvas", 800, 600)
    
    # Draw the rotated histogram
    rotated_hist.Draw("COLZ")
    
    # Calculate and draw the ellipse
    calculate_ellipse_parameters(rotated_hist, 45)

    # Update the canvas
    canvas.Draw()

# Call the function to display the rotated histogram with the ellipse
RotatedHist()