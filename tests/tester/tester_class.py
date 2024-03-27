from MDAnalysis.analysis.base import AnalysisBase


class Tester(AnalysisBase):
    """Test class for mdacli. I AM STUPID
    
    Parameters
    ----------
    bins : int
        Determines the number of bins used for data averaging; (this parameter sets the
        upper limit). The data are by default binned logarithmically. This helps to
        reduce noise, particularly in the high-frequency domain, and also prevents plot
        files from being too large.

    """

    def __init__(self, universe, **kwargs):
        """Initialise the Tester class."""
        super().__init__(universe, **kwargs)

    def _prepare(self):
        """Prepare the analysis."""
        pass

    def _single_frame(self):
        """Analyse a single frame."""
        pass

    def _conclude(self):
        """Conclude the analysis."""
        pass
