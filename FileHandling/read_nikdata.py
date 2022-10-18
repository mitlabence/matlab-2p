t301 = "D:/PhD/Data/T301/d1/anfall/T301_tmev_d1.270820.1151_nik.txt"
t370 = "D:/PhD/Data/ChR2/example with lfp/T370_ChR2_d29_elec_002_nik.txt"

t324 = "D:/PhD/Data/T324/d1/T324_mock_d1.221020.1634_nik.txt"  # imaging, untouched, utf-16 encoding.
t386 = "D:/PhD/Data/ChR2/Example nik/T386.021221.1111_nik.txt"  # stim experiment, untouched, utf-16 encoding.


def read_nikdata(fname: str = t301, encoding: str = 'utf-16'):
    """
    Open a Nikon metadata (time stamps file, mostly ending with _nik.txt) and bring it to a standardized form.
    Standardized form is how the Nikon Elements exports normal imaging data (this works well with Matlab). When using
    stimulating lasers, for example, new columns appear, and weird changes (decimal comma instead of decimal point, but
    only in some places, weird time format, new columns)
    Steps to do:
        1.  Check number of columns. Leave the file unchanged if the column names are:
                Time [s], SW Time [s], NIDAQ Time [s], Index
            If columns are:
                Time [m:s.ms], SW Time [s], NIDAQ Time [s], Index, Events, Events Type
            Then start correction.
        2.  Load data as Header and Data.
        3.  Remove rows of stimulation start and end frames. These do not correspond to imaging frames! (optional: Note
            the start and end times.)
        4.  Change Time [m:s.ms] column to Time [s] with decimal dot.
        5.  Change all commas to dots. This changes e.g. 0,1234 to 0.1234, which matlab recognises as a number.
    """
    #TODO: not all of the stimulating recordings have the weird formats!!!
    with open(fname, "r", encoding=encoding) as f:
        title = f.readline().strip()
        print(title.split("\t"))
        frow = f.readline().strip()
        print(frow.split("\t"))


read_nikdata(t370, encoding='utf-8')
