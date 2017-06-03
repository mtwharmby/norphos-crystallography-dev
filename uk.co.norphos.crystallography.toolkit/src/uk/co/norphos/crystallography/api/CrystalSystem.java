package uk.co.norphos.crystallography.api;

/**
 * An enum containing the names of the different crystal or lattice systems. 
 * 
 * The value of the enum depends on whether a crystal system or a lattice 
 * system is being described. The table below summarises which labels should 
 * be used for which case.
 * <table summary="Value to use with Lattice and Crystal System">
 * 	<tr><th>CrystalSystem (enum)	</th><th>Lattice System	</th><th>Crystal System	</th></tr>
 * 	<tr><td>CUBIC					</td><td>Cubic			</td><td>Cubic			</td></tr>
 * 	<tr><td>HEXAGONAL				</td><td>Hexagonal		</td><td>Hexagonal		</td></tr>
 * 	<tr><td>RHOMBOHEDRAL			</td><td>Rhombohedral	</td><td>Trigonal		</td></tr>
 * 	<tr><td>TRIGONAL				</td><td>Hexagonal		</td><td>Trigonal		</td></tr>
 * 	<tr><td>TETRAGONAL				</td><td>Tetragonal		</td><td>Tetragonal		</td></tr>
 * 	<tr><td>ORTHORHOMBIC			</td><td>Orthorhombic	</td><td>Orthrhombic	</td></tr>
 * 	<tr><td>MONOCLINIC				</td><td>Monoclinic		</td><td>Monoclinic		</td></tr>
 * 	<tr><td>TRICLINIC				</td><td>Triclinic		</td><td>Triclinic		</td></tr>
 * </table
 * 
 * @author Michael Wharmby
 *
 */
public enum CrystalSystem {
	TRICLINIC, MONOCLINIC, ORTHORHOMBIC, TETRAGONAL,
	TRIGONAL, RHOMBOHEDRAL, HEXAGONAL, CUBIC; 

}
