import pandas as pd
import sys
from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.live import Live

from .utilities import * 


class ProgressTable:
    """
    A class to manage and display a progress table for BlooMine runs using Rich.
    """


    def __init__(self, sample_ids, flank_ids):
        """
        Initialize the progress table.

        Args:
            sample_ids (list): A list of sample IDs (fastq prefixes).
            flank_ids (list): A list of flank IDs.
        """
        ## Create a multi-index for the table
        self.table = pd.DataFrame(index=sample_ids, columns=pd.MultiIndex.from_product([flank_ids, ['Flank 1', 'Flank 2', 'MOI']]))
        self.table.fillna(" ", inplace=True)
        self.clearlines = ((len(sample_ids)+1) * (len(flank_ids)+1)) +2
        self.console = Console()  ## Initialize Rich console
        self.printed = False  ## Flag to track if the table has been printed
        self.live = None  ## Initialize Live object for dynamic updates


    def update(self, sample_id, flank_id, step, value, color="yellow"):
        """
        Update the progress table with a new value and color.

        Args:
            sample_id (str): The sample ID.
            flank_id (str): The flank ID.
            step (str): The step being updated ('F1', 'F2', 'MOI').
            value (int or str): The new value for the step.
            color (str): The color to use for the cell ('red', 'green', 'yellow'). Default is 'yellow'.
        """
        ## Apply color based on the provided color argument
        if color == "red":
            value = f"[red]{value}[/]"  ## Red
        elif color == "green":
            value = f"[green]{value}[/]"  ## Green
        elif color == "yellow":
            value = f"[yellow]{value}[/]"  ## Yellow

        self.table.loc[sample_id, (flank_id, step)] = value
        self.print_table()


    def print_table(self, color="yellow"):
        """
        Print the current progress table to the console using Rich.
        """

        table = Table(title="BlooMine Progress")
        table.add_column("Sample ID", style="cyan", no_wrap=True)
        table.add_column("Flank ID", style="cyan", no_wrap=True)
        table.add_column("Flank 1", style="cyan", no_wrap=True)
        table.add_column("Flank 2", style="cyan", no_wrap=True)
        table.add_column("MOI", style="cyan", no_wrap=True)

        for sample_id in self.table.index:
            for flank_id in self.table.columns.levels[0]:
                table.add_row(
                    sample_id,
                    flank_id,
                    self.table.loc[sample_id, (flank_id, "Flank 1")],
                    self.table.loc[sample_id, (flank_id, "Flank 2")],
                    self.table.loc[sample_id, (flank_id, "MOI")],
                )

        if self.live is None:
            ## If Live object is not initialized, create it
            self.live = Live(table, console=self.console, refresh_per_second=4, vertical_overflow="visible")
            self.live.start()
        else:
            ## Update the Live object with the new table
            self.live.update(table)

        self.printed = True

    def close(self):
        """
        Close the Live object to stop dynamic updates.
        """
        if self.live is not None:
            self.live.stop()


def create_group(progress, table):
    progress_panel = Panel(
        progress,
        title="Progress",
        border_style="green",
        padding=(1, 1),
        height=5,
    )
    return Group(table, progress_panel)