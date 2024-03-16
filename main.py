import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import CubicSpline


class main():
    def __init__(self):
        self.open_index = 0
        self.surge_index = 4
        self.manometer_angle = 78
        self.Datum = 6.6
        self.Venturi_inlet = 0.164
        self.Venturi_outlet = 0.0507
        self.Density = 1.21559
        self.T1 = 290
        self.Surge_Margin = [0.1,0.2]
        self.Voltage = np.array([1.608679,1.637058,1.672924,1.71205,1.751427,1.790347,1.830094,1.869911,1.900875,1.873427,1.837722,1.801313,1.765257,1.729414,1.693184,1.657564,1.623563,1.590641])
        self.Mass = np.array([1,2,3,4,5,6,7,8,9,8,7,6,5,4,3,2,1,0])
        self.RPM = np.array([1500,
                    2250,
                    2625,
                    3000])
        self.theta_0 = np.radians(30)
        self.theta_1 = np.radians(50)
        self.Motor_work_rate = np.array([[1855.11724],
                                         [4510.67600],
                                         [7083.79406],
                                         [10203.1002]])
        self.Raw_Data = np.array([[
                                   [9.5,7.5,7.25,6.3,  6, 5.3,  5,  5, 9.9],
                                   [9.2,6.9,   6,4.4,3.6, 3.2,1.4,1.1, 5.2],
                                   [9.4,  7, 6.4,  5,  5, 3.2,2.4,2.4,   7],
                                   [9.5,  0,   0,  0,  0,   0,3.6,3.8, 8.4],
                                   [  9,5.8, 4.5,  2,0.5,-0.5, -3, -3,-0.2]
                                  ],
                                  [
                                   [  14,10.2,9.8, 7.8, 6.8,    5,    4,  4.4,14.8],
                                   [14.2,   0,  0,   0,   0,    0, -0.2, -0.6, 9.6],
                                   [14.2,   0,  0,   0,   0,    0,  1.6,  1.2,12.2],
                                   [14.2,   0,  0,   0,   0,    0,  2.4,    2,13.3],
                                   [12.8, 5.8,  3,-2.4,-5.2,-10.4,-13.2,-13.2,-7.4]
                                  ],
                                  [
                                   [17.2,12.3,11.2,8.6,7.4,5,3.6,3.6,18.1],
                                   [16.6,0,0,0,0,0,-13,-12.3,1],
                                   [17.2,0,0,0,0,0,-6,-6.8,6.2],
                                   [17.4,0,0,0,0,0,-2,-2.4,11.8],
                                   [15,5.8,1.6,-6.2,-9.8,-16.8,-20.8,-21,-13]
                                  ],
                                  [
                                   [21.1,14.7,13.7,9.8,8.1,4.7,3,1.5,22.5],
                                   [22.5,0,0,0,0,0,-14,-14.5,1.8],
                                   [20.5,0,0,0,0,0,-7.9,-8.5,8.7],
                                   [21.5,0,0,0,0,0,-1.5,-2.1,17.3],
                                   [18.5,6,0.8,-10,-15,-24,-27.5,-29.9,-19]
                                  ]
                                 ])
        self.cmap = plt.cm.get_cmap('turbo')
        self.colours = [self.cmap(i) for i in np.linspace(0,1,4)]
        #self.Voltage_Calibration()
        self.mainloop()
    def Voltage_Calibration(self):
        plt.scatter(self.Mass, self.Voltage)
        self.Coefficients, _ = curve_fit(self.curve_model, self.Mass, self.Voltage, p0 = [1, 1])
        self.a, self.b = self.Coefficients
        self.fitline = np.array([[x, self.curve_model(x, self.a, self.b)] for x in np.linspace(0, 10)])
        plt.plot(self.fitline[:,0],self.fitline[:,1], color = "red", label = str(self.b)+"x+"+str(self.a)+" Fit line")
        plt.legend()
        plt.xlabel("Calibration Mass [Kg]")
        plt.ylabel("Measured Voltage [Volts]")
        plt.title("Load Cell Calibration")
        plt.show()
    def curve_model(self,x, a, b):
        '''Defines the model used during the curve fitting step'''
        return (a*x)+b
    def fit_line(self, array, lower_extrapolation_factor = 1, upper_extrapolation_factor = 1):
        print(array)
        '''Find cartesian location of a point a percentage of a way along a foils chamber line using NLSR'''
        self.Coefficients, _ = curve_fit(self.curve_model, array[:,0], array[:,1], p0 = [1, 1]) #Computes the optimum Coefficients to fit the curve model to the foils points using NLSS
        print(self.Coefficients)
        self.a, self.b = self.Coefficients
        self.Sample_points = np.array([[x/100000, self.curve_model(x/100000,self.a, self.b)] for x in range(int(min(array[:,0])*100000*lower_extrapolation_factor),int(max(array[:,0])*100000*upper_extrapolation_factor))])
        return self.Sample_points
    def Work_rate(self, massflow):
        return np.min(self.Motor_work_rate)+(np.max(self.Motor_work_rate)-np.min(self.Motor_work_rate))*((massflow-np.min(self.equivalent_mass_flow_rate))/(np.max(self.equivalent_mass_flow_rate)-np.min(self.equivalent_mass_flow_rate)))

    def Free_Vortex_Design(self,i):
        self.CS = CubicSpline(self.RPM,self.Pressure_data[:,:,0][:,i])
        self.massflow_array = []
        self.compression_array = []
        for rpm in range(1,4000):
            self.U = (rpm/60)*2*np.pi*(0.420+0.300)*0.5
            self.inlet_pressure = self.CS(rpm)
            self.compressor_intake_velocity = np.sqrt(2*(-self.inlet_pressure)/self.Density)
            self.mass_flow_rate = self.compressor_intake_velocity*self.Density*0.0678
            self.phi = self.compressor_intake_velocity/self.U
            self.psi = self.phi*(np.tan(self.theta_1)-np.tan(self.theta_0))
            self.delta_h0 = self.psi*self.U**2
            self.delta_t0 = (3*self.delta_h0)/1005
            self.compression_ratio = ((self.T1+self.delta_t0)/self.T1)**(1.4/0.4)
            print("RPM: ",rpm," Temperature Change: ",self.delta_t0, "Intake Mas flow: ", self.mass_flow_rate, " Compression ratio: ", self.compression_ratio)
            self.compression_array.append(self.compression_ratio)
            self.massflow_array.append(self.mass_flow_rate)
        self.fv_Compression_ratio = np.array(self.compression_array)
        #self.fv_Compression_ratio = ((self.fv_Compression_ratio-1)/10)+1
        self.fv_massflow_array = np.array(self.massflow_array)
        plt.plot(self.fv_massflow_array, self.fv_Compression_ratio, label = "Open Flow Free Votrex Design")


    def mainloop(self):
        self.Corrected_data = (6.6-self.Raw_Data)*np.cos(np.radians(90-self.manometer_angle))
        self.Pressure_data = (self.Corrected_data/100)*9.81*784
        self.Venturi_velocity = np.sqrt((2*(self.Pressure_data[:,:,7]-self.Pressure_data[:,:,8]))/(self.Density*(1-(self.Venturi_outlet/self.Venturi_inlet)**2))   )
        self.massflow = self.Venturi_velocity*self.Venturi_outlet*self.Density
        self.temperature_change = (self.Motor_work_rate/self.massflow)/1005
        self.Compressor_efficiency = (self.massflow*(self.Pressure_data[:,:,6]-self.Pressure_data[:,:,0]))/(self.Density*self.Motor_work_rate)
        self.Compression_ratio = (self.Pressure_data[:,:,6]+100476)/(self.Pressure_data[:,:,0]+100476)
        self.equivalent_mass_flow_rate = self.massflow
        #self.Free_Vortex_Design()

        print("Raw Data")
        print(self.Raw_Data)
        print("Corrected Data")
        print(self.Corrected_data)
        print("Pressure data")
        print(self.Pressure_data)
        print("Venturi velocity")
        print(self.Venturi_velocity)
        print("Mass Flow")
        print(self.massflow)
        print("temperature change")
        print(self.temperature_change)
        print("compressor efficiency")
        print(self.Compressor_efficiency)
        print("Compression ratio")
        print(self.Compression_ratio)
        print("equivalent_mass_flow_rate")
        print(self.equivalent_mass_flow_rate)
        input()

        '''Compressor Mapping'''
        '''Plotting Raw Data'''

        '''Plotting an efficiency contour'''
        '''##############################################################################################################################################################'''
        self.intake_pressure = 100476
        self.efficiency_array = []
        self.massflow_array = np.linspace(1.14,3.4)
        self.pressure_ratio_array = np.linspace(1,1.044)
        self.massflow_array, self.pressure_ratio_array = np.meshgrid(self.massflow_array, self.pressure_ratio_array)
        self.efficiency_point =  (self.massflow_array*(self.pressure_ratio_array*100476-100476))/(self.Density*self.Work_rate(self.massflow_array))
        plt.contourf(self.massflow_array, self.pressure_ratio_array, self.efficiency_point,10, cmap='gist_yarg')
        self.colour_bar = plt.colorbar()
        self.colour_bar.set_label('Compressor Efficiency', rotation=270)




        '''Raw Data'''
        '''##############################################################################################################################################################'''
        for i in range(0,4):
            plt.scatter(self.equivalent_mass_flow_rate[i,:],self.Compression_ratio[i,:], label = str(self.RPM[i])+" RPM", color = self.colours[i])

        '''Plotting the surge line'''
        '''##############################################################################################################################################################'''
        self.Surge_line = self.fit_line(np.column_stack((self.equivalent_mass_flow_rate[:,4],self.Compression_ratio[:,4])),0.95,1.5)
        plt.plot(self.Surge_line[:,0], self.Surge_line[:,1], label = "Surge Line", color = "Red")
        plt.fill_between(self.Surge_line[:,0],self.Surge_line[:,1], 2,facecolor="none",hatch="X", label = "Compressor Stall", edgecolor="black")

        '''Working Line'''
        plt.plot(self.Surge_line[:,0], self.Surge_line[:,1]-(self.Surge_Margin[1]*(self.Surge_line[:,1]-np.min(self.Compression_ratio))), label = "20% SM Working Line", linestyle = '--')

        '''Plotting the RPM best fit lines'''
        '''##############################################################################################################################################################'''
        for i in range(0,4):
            self.line = self.fit_line(np.column_stack((self.equivalent_mass_flow_rate[i,:],self.Compression_ratio[i,:])))
            plt.plot(self.line[:,0], self.line[:,1], label = str(self.RPM[i])+" RPM fit line",color = self.colours[i])

        '''Free Vortex Line'''
        #for i in range(0,5):
        self.Free_Vortex_Design(0)


        '''Formating'''
        plt.legend()
        plt.xlabel("Equivalent Mass Flow $KgS^{-1}$")
        plt.ylabel("Compression Ratio P1/P0")
        plt.xlim(1.14, 3.4)
        plt.ylim(1,1.044)
        plt.title("Compressor Map")
        plt.show()



        plt.plot(self.massflow, self.temperature_change)
        plt.xlabel("Mass Flow")
        plt.ylabel("Temperature change")
        plt.show()

        plt.plot(self.fv_massflow_array, self.fv_Compression_ratio)
        plt.show()




if __name__ == '__main__':
    main()
