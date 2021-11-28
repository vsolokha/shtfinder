from kivy.app import App
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.uix.togglebutton import ToggleButton
import shtfinder as sf
import numpy as np
import webbrowser


class GUIScreen(GridLayout):

    def __init__(self, **kwargs):
        super(GUIScreen, self).__init__(**kwargs)
        self.cols = 2
        self.unpack_method = "exe"
        self.add_widget(Label(text='BASESHOT'))
        self.baseshot = TextInput(multiline=False, text='36612')
        self.add_widget(self.baseshot)
        self.add_widget(Label(text='START SEQUENCE'))
        self.start_sequence = TextInput(multiline=False, text='36610')
        self.add_widget(self.start_sequence)
        self.add_widget(Label(text='STOP SEQUENCE'))
        self.stop_sequence = TextInput(multiline=False, text='36614')
        self.add_widget(self.stop_sequence)
        self.add_widget(Label(text='DIAGNOSTICS'))
        self.diag = TextInput(multiline=False, text='6,7')
        self.add_widget(self.diag)
        self.add_widget(Label(text='WEIGHTS'))
        self.weights = TextInput(multiline=False, text='1.0,0.0')
        self.add_widget(self.weights)
        self.add_widget(Label(text='SHT PATH'))
        self.shtpath = TextInput(multiline=False, text='./data')
        self.add_widget(self.shtpath)
        self.btn1 = Button(text='Show Diagnostic Names')
        self.add_widget(self.btn1)
        self.btn2 = Button(text='Get Distances')
        self.add_widget(self.btn2)
        self.btn3 = ToggleButton(text='Use ShtRipper (slow)', state='normal')
        self.add_widget(self.btn3)
        self.diag_info = Label(text='DIAGNOSTICS \n INFO')
        self.add_widget(self.diag_info)
        self.set_buttons()
        
    def set_buttons(self):
        self.btn1.bind(on_press=self.show_diag)
        self.btn2.bind(on_press=self.get_distances)
        self.btn3.bind(on_press=self.choose_ripper_unpack)

    def get_distances(self, instance):
        self.args = [self.baseshot.text, self.start_sequence.text, self.stop_sequence.text, self.diag.text, self.weights.text, self.shtpath.text]
        baseshot   = int(self.args[0])
        test_range = [int(self.args[1]), int(self.args[2])] 
        columns    = list(map(int,   self.args[3].split(",")))
        weights    = list(map(float, self.args[4].split(",")))
        shtpath    = str(self.args[5])
        res_dtw, __     = sf.get_distance(baseshot, test_range, columns, shtpath, "dtw", self.unpack_method)
        res_euc, number = sf.get_distance(baseshot, test_range, columns, shtpath, "euc", self.unpack_method)
        res_dtw_sort = sf.get_wsorted(number, res_dtw, weights)
        print(res_dtw_sort)
        try:
            np.savetxt(f"output/result-dtw-{test_range[0]}.txt", res_dtw_sort, \
                fmt=['%6d','%5.4e'], header='DTW DISTANCE')
        except: 
            np.savetxt(f"output/result-dtw-{test_range[0]}.txt", res_dtw_sort, \
                header='DTW DISTANCE') 
        res_euc_sort = sf.get_wsorted(number, res_euc, weights)
        try:
            np.savetxt(f"output/result-euc-{test_range[0]}.txt", res_euc_sort, \
                fmt=['%6d','%5.4e'], header='EUCLIDIAN DISTANCE')
        except:
            np.savetxt(f"output/result-euc-{test_range[0]}.txt", res_euc_sort, \
                header='EUCLIDIAN DISTANCE')
        webbrowser.open(f"output/result-dtw-{test_range[0]}.txt")
        webbrowser.open(f"output/result-euc-{test_range[0]}.txt")    

    def show_diag(self, instance):
        self.args = [self.baseshot.text, self.start_sequence.text, self.stop_sequence.text, self.diag.text, self.weights.text, self.shtpath.text]
        baseshot_num = int(self.args[0])
        shtpath  = str(self.args[5])
        baseshot = sf.Shot(baseshot_num, shtpath, self.unpack_method) 
        i = 0
        names = ""
        for name in baseshot.names:
            names += str(i) + "\t" + name + "\n"
            i += 1
        with open("./output/names.txt", "w") as tf:
            tf.write(names)
        webbrowser.open("output/names.txt")    
        self.diag_info.text = "Names are in the WebBrowser Window"
        
    def choose_ripper_unpack(self, instance):
        if instance.state == "down":
            self.unpack_method = "shtripper"
        else:
            self.unpack_method = "exe"

       
class MyApp(App):

    def build(self):
        return GUIScreen()


if __name__ == '__main__':
    MyApp().run()
