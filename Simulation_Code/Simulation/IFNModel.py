from cc3d import CompuCellSetup
        
from IFNModelSteppables import ODEModelSteppable
CompuCellSetup.register_steppable(steppable=ODEModelSteppable(frequency=1))

from IFNModelSteppables import CellularModelSteppable
CompuCellSetup.register_steppable(steppable=CellularModelSteppable(frequency=1))

from IFNModelSteppables import PlaqueAssaySteppable
CompuCellSetup.register_steppable(steppable=PlaqueAssaySteppable(frequency=1))

from IFNModelSteppables import OutputSteppable
CompuCellSetup.register_steppable(steppable=OutputSteppable(frequency=1))

CompuCellSetup.run()
