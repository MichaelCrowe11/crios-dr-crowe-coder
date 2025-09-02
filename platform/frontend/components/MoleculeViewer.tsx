// 3D Molecule Viewer Component
// Interactive Three.js molecular visualization

import React, { useRef, useEffect, useState } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Text, Sphere, Cylinder } from '@react-three/drei';
import * as THREE from 'three';
import { Box } from '@mui/material';

interface Atom {
  id: string;
  element: string;
  position: [number, number, number];
}

interface Bond {
  from: string;
  to: string;
  order: number;
}

interface MoleculeViewerProps {
  smiles: string;
}

// Atom colors based on CPK coloring
const atomColors: Record<string, string> = {
  H: '#ffffff',
  C: '#909090',
  N: '#3050f8',
  O: '#ff0d0d',
  F: '#90e050',
  P: '#ff8000',
  S: '#ffff30',
  Cl: '#1ff01f',
  Br: '#a62929',
  I: '#940094'
};

// Atom sizes (van der Waals radii)
const atomSizes: Record<string, number> = {
  H: 0.25,
  C: 0.35,
  N: 0.32,
  O: 0.30,
  F: 0.28,
  P: 0.38,
  S: 0.38,
  Cl: 0.35,
  Br: 0.38,
  I: 0.43
};

function MoleculeModel({ atoms, bonds }: { atoms: Atom[], bonds: Bond[] }) {
  const groupRef = useRef<THREE.Group>(null);

  useFrame((state, delta) => {
    if (groupRef.current) {
      groupRef.current.rotation.y += delta * 0.2;
    }
  });

  const getAtomPosition = (atomId: string): [number, number, number] => {
    const atom = atoms.find(a => a.id === atomId);
    return atom ? atom.position : [0, 0, 0];
  };

  const getBondPosition = (bond: Bond): [number, number, number] => {
    const pos1 = getAtomPosition(bond.from);
    const pos2 = getAtomPosition(bond.to);
    return [
      (pos1[0] + pos2[0]) / 2,
      (pos1[1] + pos2[1]) / 2,
      (pos1[2] + pos2[2]) / 2
    ];
  };

  const getBondLength = (bond: Bond): number => {
    const pos1 = getAtomPosition(bond.from);
    const pos2 = getAtomPosition(bond.to);
    const dx = pos2[0] - pos1[0];
    const dy = pos2[1] - pos1[1];
    const dz = pos2[2] - pos1[2];
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
  };

  const getBondRotation = (bond: Bond): THREE.Euler => {
    const pos1 = getAtomPosition(bond.from);
    const pos2 = getAtomPosition(bond.to);
    const direction = new THREE.Vector3(
      pos2[0] - pos1[0],
      pos2[1] - pos1[1],
      pos2[2] - pos1[2]
    ).normalize();
    
    const quaternion = new THREE.Quaternion();
    quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction);
    const euler = new THREE.Euler().setFromQuaternion(quaternion);
    return euler;
  };

  return (
    <group ref={groupRef}>
      {/* Render atoms */}
      {atoms.map((atom) => (
        <Sphere
          key={atom.id}
          position={atom.position}
          args={[atomSizes[atom.element] || 0.3, 32, 32]}
        >
          <meshPhongMaterial 
            color={atomColors[atom.element] || '#808080'} 
            emissive={atomColors[atom.element] || '#808080'}
            emissiveIntensity={0.2}
            shininess={100}
          />
        </Sphere>
      ))}

      {/* Render bonds */}
      {bonds.map((bond, index) => {
        const position = getBondPosition(bond);
        const length = getBondLength(bond);
        const rotation = getBondRotation(bond);

        if (bond.order === 1) {
          return (
            <Cylinder
              key={`bond-${index}`}
              position={position}
              rotation={rotation}
              args={[0.1, 0.1, length, 8]}
            >
              <meshPhongMaterial color="#666666" />
            </Cylinder>
          );
        } else if (bond.order === 2) {
          return (
            <group key={`bond-${index}`}>
              <Cylinder
                position={[position[0] - 0.1, position[1], position[2]]}
                rotation={rotation}
                args={[0.08, 0.08, length, 8]}
              >
                <meshPhongMaterial color="#666666" />
              </Cylinder>
              <Cylinder
                position={[position[0] + 0.1, position[1], position[2]]}
                rotation={rotation}
                args={[0.08, 0.08, length, 8]}
              >
                <meshPhongMaterial color="#666666" />
              </Cylinder>
            </group>
          );
        } else {
          return (
            <group key={`bond-${index}`}>
              <Cylinder
                position={[position[0] - 0.15, position[1], position[2]]}
                rotation={rotation}
                args={[0.06, 0.06, length, 8]}
              >
                <meshPhongMaterial color="#666666" />
              </Cylinder>
              <Cylinder
                position={position}
                rotation={rotation}
                args={[0.06, 0.06, length, 8]}
              >
                <meshPhongMaterial color="#666666" />
              </Cylinder>
              <Cylinder
                position={[position[0] + 0.15, position[1], position[2]]}
                rotation={rotation}
                args={[0.06, 0.06, length, 8]}
              >
                <meshPhongMaterial color="#666666" />
              </Cylinder>
            </group>
          );
        }
      })}
    </group>
  );
}

const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ smiles }) => {
  const [atoms, setAtoms] = useState<Atom[]>([]);
  const [bonds, setBonds] = useState<Bond[]>([]);

  useEffect(() => {
    // Parse SMILES and generate 3D coordinates
    // For demo, using aspirin structure
    if (smiles === 'CC(=O)Oc1ccccc1C(=O)O') {
      // Aspirin structure
      setAtoms([
        { id: 'C1', element: 'C', position: [0, 0, 0] },
        { id: 'C2', element: 'C', position: [1.5, 0, 0] },
        { id: 'O1', element: 'O', position: [2.5, 1, 0] },
        { id: 'O2', element: 'O', position: [2.5, -1, 0] },
        { id: 'C3', element: 'C', position: [-1.5, 0, 0] },
        { id: 'C4', element: 'C', position: [-2.5, 1, 0] },
        { id: 'C5', element: 'C', position: [-3.8, 1, 0] },
        { id: 'C6', element: 'C', position: [-4.5, 0, 0] },
        { id: 'C7', element: 'C', position: [-3.8, -1, 0] },
        { id: 'C8', element: 'C', position: [-2.5, -1, 0] },
        { id: 'C9', element: 'C', position: [-2.5, -2.5, 0] },
        { id: 'O3', element: 'O', position: [-1.5, -3.5, 0] },
        { id: 'O4', element: 'O', position: [-3.5, -2.5, 0] }
      ]);

      setBonds([
        { from: 'C1', to: 'C2', order: 1 },
        { from: 'C2', to: 'O1', order: 2 },
        { from: 'C2', to: 'O2', order: 1 },
        { from: 'O2', to: 'C3', order: 1 },
        { from: 'C3', to: 'C4', order: 2 },
        { from: 'C4', to: 'C5', order: 1 },
        { from: 'C5', to: 'C6', order: 2 },
        { from: 'C6', to: 'C7', order: 1 },
        { from: 'C7', to: 'C8', order: 2 },
        { from: 'C8', to: 'C3', order: 1 },
        { from: 'C8', to: 'C9', order: 1 },
        { from: 'C9', to: 'O3', order: 2 },
        { from: 'C9', to: 'O4', order: 1 }
      ]);
    }
  }, [smiles]);

  return (
    <Box sx={{ width: '100%', height: '100%', bgcolor: '#000' }}>
      <Canvas camera={{ position: [0, 0, 10], fov: 60 }}>
        <ambientLight intensity={0.5} />
        <pointLight position={[10, 10, 10]} intensity={1} />
        <pointLight position={[-10, -10, -10]} intensity={0.5} />
        <MoleculeModel atoms={atoms} bonds={bonds} />
        <OrbitControls enablePan={true} enableZoom={true} enableRotate={true} />
      </Canvas>
    </Box>
  );
};

export default MoleculeViewer;