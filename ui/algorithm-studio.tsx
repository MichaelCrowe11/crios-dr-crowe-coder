import React, { useState, useEffect, useRef, useCallback } from "react";
import { motion, AnimatePresence } from "framer-motion";
import {
  Card,
  CardContent,
  CardHeader,
  CardTitle,
} from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Input } from "@/components/ui/input";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Slider } from "@/components/ui/slider";
import { Switch } from "@/components/ui/switch";
import { Progress } from "@/components/ui/progress";
import { Textarea } from "@/components/ui/textarea";
import { ScrollArea } from "@/components/ui/scroll-area";
import { Separator } from "@/components/ui/separator";
import { Alert, AlertDescription } from "@/components/ui/alert";
import { 
  Code2, 
  Sparkles, 
  Cpu, 
  Dna, 
  FlaskConical, 
  Brain, 
  Zap, 
  Play, 
  Pause, 
  Save, 
  Download, 
  Upload,
  GitBranch,
  Layers,
  Activity,
  Microscope,
  Atom,
  Binary,
  Network,
  Shield,
  Gauge,
  Workflow,
  TestTube,
  LineChart,
  Bot,
  Lightbulb,
  Rocket,
  Eye,
  Copy,
  Check,
  X,
  AlertCircle,
  Info,
  ChevronRight,
  ChevronDown,
  Plus,
  Minus,
  RefreshCw,
  Settings,
  Terminal,
  FileCode,
  Database,
  Cloud,
  Lock,
  Unlock,
  Timer,
  TrendingUp,
  Target,
  Beaker,
  Shuffle,
  GitMerge,
  Package,
  Boxes,
  Fingerprint
} from "lucide-react";
import Editor from "@monaco-editor/react";
import { ResponsiveContainer, LineChart, Line, AreaChart, Area, XAxis, YAxis, CartesianGrid, Tooltip, Legend, RadarChart, PolarGrid, PolarAngleAxis, PolarRadiusAxis, Radar } from "recharts";

// ------------------------------
// Algorithm Generation Engine
// ------------------------------

interface AlgorithmBlueprint {
  id: string;
  name: string;
  category: "discovery" | "optimization" | "validation" | "synthesis" | "safety";
  description: string;
  inputs: AlgorithmInput[];
  outputs: AlgorithmOutput[];
  code: string;
  performance: PerformanceMetrics;
  biologicalPatterns: string[];
  ethicsScore: number;
}

interface AlgorithmInput {
  name: string;
  type: "molecule" | "protein" | "dataset" | "parameters";
  format: string;
  required: boolean;
}

interface AlgorithmOutput {
  name: string;
  type: "compounds" | "scores" | "predictions" | "visualizations";
  format: string;
}

interface PerformanceMetrics {
  speed: number; // ms
  accuracy: number; // 0-1
  scalability: number; // 0-1
  memory: number; // MB
}

interface GenerationRequest {
  objective: string;
  constraints: string[];
  targetMetrics: {
    minAccuracy: number;
    maxLatency: number;
    ethicsCompliance: boolean;
  };
  biologicalInspiration: boolean;
  quantumReady: boolean;
}

// Visual Node Editor Types
interface AlgorithmNode {
  id: string;
  type: "input" | "process" | "output" | "condition" | "loop";
  label: string;
  position: { x: number; y: number };
  data: any;
  connections: string[];
}

interface NodeConnection {
  from: string;
  to: string;
  type: "data" | "control";
}

// AI Agent for Algorithm Generation
class AlgorithmGenerator {
  private templates: Map<string, string>;
  private patterns: Map<string, string>;
  
  constructor() {
    this.templates = new Map();
    this.patterns = new Map();
    this.initializeTemplates();
    this.initializePatterns();
  }

  private initializeTemplates() {
    // Drug Discovery Templates
    this.templates.set("lead_optimization", `
def optimize_lead_compound(smiles: str, target: str) -> List[Compound]:
    """
    CriOS Lead Optimization Algorithm
    Applies Crowe Framework for multi-objective optimization
    """
    mol = Chem.MolFromSmiles(smiles)
    
    # Generate analogs using bioisosteric replacements
    analogs = []
    replacements = get_bioisosteric_replacements(mol)
    
    for replacement in replacements:
        analog = apply_replacement(mol, replacement)
        
        # Calculate Crowe Score components
        potency = predict_potency(analog, target)
        selectivity = predict_selectivity(analog, target)
        admet = predict_admet_properties(analog)
        synthesis = assess_synthetic_accessibility(analog)
        safety = predict_safety_profile(analog)
        
        # Apply ethics filter
        if passes_ethics_screen(analog):
            score = calculate_crowe_score(
                potency, selectivity, admet, synthesis, safety
            )
            analogs.append({
                'smiles': Chem.MolToSmiles(analog),
                'crowe_score': score,
                'properties': extract_properties(analog)
            })
    
    return sorted(analogs, key=lambda x: x['crowe_score'], reverse=True)[:20]
    `);

    this.templates.set("target_prediction", `
def predict_protein_targets(smiles: str, threshold: float = 0.7) -> List[Target]:
    """
    AI-Powered Target Prediction using Graph Neural Networks
    """
    mol = Chem.MolFromSmiles(smiles)
    
    # Generate molecular graph representation
    graph = molecule_to_graph(mol)
    
    # Load pre-trained GNN model
    model = load_model('crios_target_predictor_v2.pkl')
    
    # Predict binding probabilities
    predictions = model.predict(graph)
    
    # Filter by threshold and ethics
    targets = []
    for protein, probability in predictions.items():
        if probability >= threshold:
            if is_ethical_target(protein):
                targets.append({
                    'protein': protein,
                    'probability': probability,
                    'mechanism': predict_mechanism(mol, protein),
                    'safety_flags': check_off_targets(mol, protein)
                })
    
    return sorted(targets, key=lambda x: x['probability'], reverse=True)
    `);

    this.templates.set("quantum_docking", `
def quantum_enhanced_docking(ligand: str, protein: str) -> DockingResult:
    """
    Quantum-Classical Hybrid Docking Algorithm
    Leverages quantum annealing for conformational search
    """
    # Classical preprocessing
    ligand_mol = prepare_ligand(ligand)
    protein_structure = load_protein(protein)
    binding_site = identify_binding_site(protein_structure)
    
    # Quantum conformational sampling
    conformers = quantum_conformer_generation(ligand_mol, n_qubits=16)
    
    # Classical scoring with ML enhancement
    docking_poses = []
    for conformer in conformers:
        pose = dock_conformer(conformer, binding_site)
        
        # Hybrid scoring function
        classical_score = calculate_glide_score(pose)
        ml_score = predict_binding_affinity(pose)
        quantum_correction = quantum_interaction_energy(pose)
        
        total_score = weighted_average([
            classical_score, ml_score, quantum_correction
        ])
        
        docking_poses.append({
            'pose': pose,
            'score': total_score,
            'interactions': analyze_interactions(pose)
        })
    
    return select_best_pose(docking_poses)
    `);
  }

  private initializePatterns() {
    this.patterns.set("homeostasis", "self.regulate_parameters()");
    this.patterns.set("evolution", "population = evolve_population(population, fitness_fn)");
    this.patterns.set("swarm", "solutions = particle_swarm_optimization(objective, bounds)");
    this.patterns.set("fractal", "structure = generate_fractal_scaffold(seed, iterations)");
    this.patterns.set("neural", "output = neural_network.forward(input)");
  }

  generateAlgorithm(request: GenerationRequest): string {
    // AI-powered algorithm generation
    let code = this.selectTemplate(request.objective);
    
    if (request.biologicalInspiration) {
      code = this.applyBiologicalPatterns(code);
    }
    
    if (request.quantumReady) {
      code = this.addQuantumEnhancements(code);
    }
    
    code = this.enforceEthicsConstraints(code, request.constraints);
    code = this.optimizePerformance(code, request.targetMetrics);
    
    return code;
  }

  private selectTemplate(objective: string): string {
    // Smart template selection based on objective
    if (objective.includes("optimization")) {
      return this.templates.get("lead_optimization") || "";
    } else if (objective.includes("target")) {
      return this.templates.get("target_prediction") || "";
    } else if (objective.includes("docking")) {
      return this.templates.get("quantum_docking") || "";
    }
    return "";
  }

  private applyBiologicalPatterns(code: string): string {
    // Inject biological computing patterns
    return code.replace(
      /optimize\(/g,
      "bio_optimize("
    ).replace(
      /predict\(/g, 
      "neural_predict("
    );
  }

  private addQuantumEnhancements(code: string): string {
    // Add quantum computing capabilities
    return `from qiskit import QuantumCircuit, execute\n${code}`;
  }

  private enforceEthicsConstraints(code: string, constraints: string[]): string {
    // Add ethics checks
    const ethicsCheck = `
    # CriOS Ethics Enforcement
    if not passes_ethics_screen(compound):
        continue  # Skip unethical compounds
    `;
    return code.replace(/return /g, `${ethicsCheck}\n    return `);
  }

  private optimizePerformance(code: string, metrics: any): string {
    // Performance optimization
    if (metrics.maxLatency < 100) {
      code = `@jit\n${code}`; // Add JIT compilation
    }
    return code;
  }
}

// Visual Algorithm Builder Component
function VisualAlgorithmBuilder({ onGenerate }: { onGenerate: (code: string) => void }) {
  const [nodes, setNodes] = useState<AlgorithmNode[]>([
    {
      id: "input1",
      type: "input",
      label: "SMILES Input",
      position: { x: 50, y: 150 },
      data: { format: "smiles" },
      connections: ["process1"]
    },
    {
      id: "process1",
      type: "process",
      label: "Molecular Processing",
      position: { x: 250, y: 150 },
      data: { algorithm: "rdkit" },
      connections: ["process2"]
    },
    {
      id: "process2",
      type: "process",
      label: "Crowe Scoring",
      position: { x: 450, y: 150 },
      data: { algorithm: "crowe_framework" },
      connections: ["output1"]
    },
    {
      id: "output1",
      type: "output",
      label: "Scored Compounds",
      position: { x: 650, y: 150 },
      data: { format: "json" },
      connections: []
    }
  ]);

  const [selectedNode, setSelectedNode] = useState<string | null>(null);
  const [dragging, setDragging] = useState<string | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);

  const handleNodeDrag = (nodeId: string, newPosition: { x: number; y: number }) => {
    setNodes(prev => prev.map(node => 
      node.id === nodeId ? { ...node, position: newPosition } : node
    ));
  };

  const addNode = (type: AlgorithmNode["type"]) => {
    const newNode: AlgorithmNode = {
      id: `node_${Date.now()}`,
      type,
      label: `New ${type}`,
      position: { x: 350, y: 250 },
      data: {},
      connections: []
    };
    setNodes(prev => [...prev, newNode]);
  };

  const generateFromVisual = () => {
    // Convert visual flow to code
    let code = "# Generated by CriOS Algorithm Studio\n\n";
    
    nodes.forEach(node => {
      switch (node.type) {
        case "input":
          code += `# Input: ${node.label}\n`;
          code += `input_data = load_${node.data.format}()\n\n`;
          break;
        case "process":
          code += `# Process: ${node.label}\n`;
          code += `result = ${node.data.algorithm}(input_data)\n\n`;
          break;
        case "output":
          code += `# Output: ${node.label}\n`;
          code += `save_${node.data.format}(result)\n\n`;
          break;
      }
    });
    
    onGenerate(code);
  };

  return (
    <Card className="rounded-2xl shadow-sm">
      <CardHeader>
        <CardTitle className="text-base flex items-center gap-2">
          <Workflow className="h-4 w-4" />
          Visual Algorithm Builder
        </CardTitle>
      </CardHeader>
      <CardContent>
        <div className="space-y-4">
          {/* Toolbar */}
          <div className="flex gap-2">
            <Button size="sm" onClick={() => addNode("input")} className="gap-1">
              <Plus className="h-3 w-3" /> Input
            </Button>
            <Button size="sm" onClick={() => addNode("process")} className="gap-1">
              <Plus className="h-3 w-3" /> Process
            </Button>
            <Button size="sm" onClick={() => addNode("condition")} className="gap-1">
              <Plus className="h-3 w-3" /> Condition
            </Button>
            <Button size="sm" onClick={() => addNode("output")} className="gap-1">
              <Plus className="h-3 w-3" /> Output
            </Button>
            <Separator orientation="vertical" />
            <Button size="sm" variant="default" onClick={generateFromVisual} className="gap-1">
              <Code2 className="h-3 w-3" /> Generate Code
            </Button>
          </div>

          {/* Canvas */}
          <div className="relative w-full h-96 bg-muted/20 rounded-xl border-2 border-dashed border-muted overflow-hidden">
            <svg ref={svgRef} className="absolute inset-0 w-full h-full">
              {/* Draw connections */}
              {nodes.map(node => 
                node.connections.map(targetId => {
                  const target = nodes.find(n => n.id === targetId);
                  if (!target) return null;
                  return (
                    <line
                      key={`${node.id}-${targetId}`}
                      x1={node.position.x + 75}
                      y1={node.position.y + 25}
                      x2={target.position.x + 75}
                      y2={target.position.y + 25}
                      stroke="#888"
                      strokeWidth="2"
                      markerEnd="url(#arrowhead)"
                    />
                  );
                })
              )}
              <defs>
                <marker id="arrowhead" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto">
                  <polygon points="0 0, 10 3, 0 6" fill="#888" />
                </marker>
              </defs>
            </svg>

            {/* Render nodes */}
            {nodes.map(node => (
              <motion.div
                key={node.id}
                className={`absolute w-36 h-12 rounded-lg border-2 flex items-center justify-center text-xs font-medium cursor-move select-none ${
                  node.type === "input" ? "bg-blue-100 border-blue-300" :
                  node.type === "process" ? "bg-green-100 border-green-300" :
                  node.type === "condition" ? "bg-yellow-100 border-yellow-300" :
                  "bg-purple-100 border-purple-300"
                } ${selectedNode === node.id ? "ring-2 ring-primary" : ""}`}
                style={{ 
                  left: node.position.x, 
                  top: node.position.y,
                  zIndex: dragging === node.id ? 1000 : 1
                }}
                drag
                dragMomentum={false}
                onDragEnd={(e, info) => {
                  handleNodeDrag(node.id, {
                    x: node.position.x + info.offset.x,
                    y: node.position.y + info.offset.y
                  });
                }}
                onClick={() => setSelectedNode(node.id)}
                onMouseDown={() => setDragging(node.id)}
                onMouseUp={() => setDragging(null)}
                whileHover={{ scale: 1.05 }}
                whileTap={{ scale: 0.95 }}
              >
                {node.type === "input" && <Upload className="h-3 w-3 mr-1" />}
                {node.type === "process" && <Cpu className="h-3 w-3 mr-1" />}
                {node.type === "condition" && <GitBranch className="h-3 w-3 mr-1" />}
                {node.type === "output" && <Download className="h-3 w-3 mr-1" />}
                {node.label}
              </motion.div>
            ))}
          </div>

          {/* Node Properties */}
          {selectedNode && (
            <Card className="rounded-xl">
              <CardHeader className="pb-3">
                <CardTitle className="text-sm">Node Properties</CardTitle>
              </CardHeader>
              <CardContent className="space-y-2">
                <Input 
                  placeholder="Node Label" 
                  value={nodes.find(n => n.id === selectedNode)?.label || ""}
                  onChange={(e) => {
                    setNodes(prev => prev.map(node => 
                      node.id === selectedNode ? { ...node, label: e.target.value } : node
                    ));
                  }}
                />
                <Select defaultValue="rdkit">
                  <SelectTrigger>
                    <SelectValue placeholder="Algorithm" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="rdkit">RDKit Processing</SelectItem>
                    <SelectItem value="crowe_framework">Crowe Framework</SelectItem>
                    <SelectItem value="ethics_filter">Ethics Filter</SelectItem>
                    <SelectItem value="ml_predictor">ML Predictor</SelectItem>
                  </SelectContent>
                </Select>
              </CardContent>
            </Card>
          )}
        </div>
      </CardContent>
    </Card>
  );
}

// Algorithm Performance Monitor
function PerformanceMonitor({ metrics }: { metrics: PerformanceMetrics }) {
  const radarData = [
    { metric: "Speed", value: metrics.speed / 10, fullMark: 100 },
    { metric: "Accuracy", value: metrics.accuracy * 100, fullMark: 100 },
    { metric: "Scalability", value: metrics.scalability * 100, fullMark: 100 },
    { metric: "Memory", value: (1000 - metrics.memory) / 10, fullMark: 100 },
  ];

  return (
    <Card className="rounded-2xl shadow-sm">
      <CardHeader>
        <CardTitle className="text-base flex items-center gap-2">
          <Gauge className="h-4 w-4" />
          Performance Metrics
        </CardTitle>
      </CardHeader>
      <CardContent>
        <ResponsiveContainer width="100%" height={200}>
          <RadarChart data={radarData}>
            <PolarGrid />
            <PolarAngleAxis dataKey="metric" />
            <PolarRadiusAxis angle={90} domain={[0, 100]} />
            <Radar name="Performance" dataKey="value" stroke="#3b82f6" fill="#3b82f6" fillOpacity={0.6} />
          </RadarChart>
        </ResponsiveContainer>
        <div className="grid grid-cols-2 gap-2 mt-4">
          <div className="text-center">
            <div className="text-2xl font-bold text-green-600">{(metrics.accuracy * 100).toFixed(1)}%</div>
            <div className="text-xs text-muted-foreground">Accuracy</div>
          </div>
          <div className="text-center">
            <div className="text-2xl font-bold text-blue-600">{metrics.speed}ms</div>
            <div className="text-xs text-muted-foreground">Latency</div>
          </div>
        </div>
      </CardContent>
    </Card>
  );
}

// Main Algorithm Studio Component
export default function AlgorithmStudio() {
  const [code, setCode] = useState<string>("");
  const [language, setLanguage] = useState<"python" | "typescript" | "rust">("python");
  const [isGenerating, setIsGenerating] = useState(false);
  const [isRunning, setIsRunning] = useState(false);
  const [output, setOutput] = useState<string>("");
  const [selectedTemplate, setSelectedTemplate] = useState<string>("");
  const [biologicalMode, setBiologicalMode] = useState(true);
  const [quantumMode, setQuantumMode] = useState(false);
  const [ethicsEnforced, setEthicsEnforced] = useState(true);
  const [generationPrompt, setGenerationPrompt] = useState("");
  const [activeTab, setActiveTab] = useState("editor");
  const [performanceMetrics, setPerformanceMetrics] = useState<PerformanceMetrics>({
    speed: 45,
    accuracy: 0.94,
    scalability: 0.87,
    memory: 256
  });

  const generator = useRef(new AlgorithmGenerator());

  // Algorithm Templates
  const templates = [
    {
      id: "lead_opt",
      name: "Lead Optimization",
      icon: <Beaker className="h-4 w-4" />,
      description: "Multi-objective compound optimization"
    },
    {
      id: "target_pred",
      name: "Target Prediction",
      icon: <Target className="h-4 w-4" />,
      description: "AI-powered protein target identification"
    },
    {
      id: "quantum_dock",
      name: "Quantum Docking",
      icon: <Atom className="h-4 w-4" />,
      description: "Quantum-enhanced molecular docking"
    },
    {
      id: "synthesis",
      name: "Retrosynthesis",
      icon: <GitMerge className="h-4 w-4" />,
      description: "AI-driven synthesis route planning"
    },
    {
      id: "admet",
      name: "ADMET Prediction",
      icon: <Activity className="h-4 w-4" />,
      description: "Pharmacokinetic property prediction"
    },
    {
      id: "toxicity",
      name: "Toxicity Screening",
      icon: <Shield className="h-4 w-4" />,
      description: "Safety and toxicity assessment"
    }
  ];

  const generateAlgorithm = async () => {
    setIsGenerating(true);
    
    // Simulate AI generation
    setTimeout(() => {
      const request: GenerationRequest = {
        objective: generationPrompt || "optimize drug-like compounds for CNS targets",
        constraints: ethicsEnforced ? ["no_controlled_substances", "no_toxic_fragments"] : [],
        targetMetrics: {
          minAccuracy: 0.9,
          maxLatency: 100,
          ethicsCompliance: ethicsEnforced
        },
        biologicalInspiration: biologicalMode,
        quantumReady: quantumMode
      };
      
      const generatedCode = generator.current.generateAlgorithm(request);
      setCode(generatedCode);
      setIsGenerating(false);
      
      // Update performance metrics
      setPerformanceMetrics({
        speed: Math.random() * 50 + 20,
        accuracy: Math.random() * 0.2 + 0.8,
        scalability: Math.random() * 0.3 + 0.7,
        memory: Math.random() * 500 + 100
      });
    }, 2000);
  };

  const runAlgorithm = async () => {
    setIsRunning(true);
    setOutput("Initializing CriOS Algorithm Engine...\n");
    
    // Simulate execution
    const steps = [
      "Loading molecular libraries...",
      "Initializing Crowe Framework...",
      "Applying biological patterns...",
      "Running ethics filters...",
      "Optimizing compounds...",
      "Calculating scores...",
      "Generating results..."
    ];
    
    for (const step of steps) {
      await new Promise(resolve => setTimeout(resolve, 500));
      setOutput(prev => prev + step + "\n");
    }
    
    setOutput(prev => prev + "\n✅ Algorithm executed successfully!\n");
    setOutput(prev => prev + "🧬 Generated 47 optimized compounds\n");
    setOutput(prev => prev + "⚡ Average Crowe Score: 0.87\n");
    setOutput(prev => prev + "🛡️ Ethics compliance: 100%\n");
    setIsRunning(false);
  };

  return (
    <div className="p-6 md:p-8 space-y-6 max-w-[1600px] mx-auto">
      {/* Header */}
      <div className="flex flex-col md:flex-row md:items-end md:justify-between gap-4">
        <div>
          <h1 className="text-2xl md:text-3xl font-semibold tracking-tight flex items-center gap-3">
            <motion.div
              className="rounded-xl bg-gradient-to-br from-blue-500 via-purple-500 to-pink-500 p-2"
              animate={{ rotate: [0, 5, -5, 0] }}
              transition={{ duration: 2, repeat: Infinity }}
            >
              <Code2 className="h-6 w-6 text-white" />
            </motion.div>
            CriOS Algorithm Studio
          </h1>
          <p className="text-sm text-muted-foreground">AI-Powered Algorithm Generation for Drug Discovery</p>
        </div>
        <div className="flex gap-2 items-center">
          <Badge variant="outline" className="gap-1">
            <Cpu className="h-3 w-3" />
            {language.toUpperCase()}
          </Badge>
          <Badge variant={biologicalMode ? "default" : "outline"} className="gap-1">
            <Dna className="h-3 w-3" />
            Biological
          </Badge>
          <Badge variant={quantumMode ? "default" : "outline"} className="gap-1">
            <Atom className="h-3 w-3" />
            Quantum
          </Badge>
          <Badge variant={ethicsEnforced ? "default" : "outline"} className="gap-1">
            <Shield className="h-3 w-3" />
            Ethics
          </Badge>
        </div>
      </div>

      {/* Generation Controls */}
      <Card className="rounded-2xl shadow-sm">
        <CardHeader>
          <CardTitle className="text-base flex items-center gap-2">
            <Sparkles className="h-4 w-4" />
            Algorithm Generation
          </CardTitle>
        </CardHeader>
        <CardContent>
          <div className="space-y-4">
            <div className="flex gap-2">
              <Input
                placeholder="Describe the algorithm you want to generate (e.g., 'optimize kinase inhibitors with BBB penetration')"
                value={generationPrompt}
                onChange={(e) => setGenerationPrompt(e.target.value)}
                className="flex-1"
              />
              <Button 
                onClick={generateAlgorithm} 
                disabled={isGenerating}
                className="gap-2"
              >
                {isGenerating ? (
                  <>
                    <RefreshCw className="h-4 w-4 animate-spin" />
                    Generating...
                  </>
                ) : (
                  <>
                    <Sparkles className="h-4 w-4" />
                    Generate
                  </>
                )}
              </Button>
            </div>

            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-2">
              {templates.map(template => (
                <Button
                  key={template.id}
                  variant={selectedTemplate === template.id ? "default" : "outline"}
                  size="sm"
                  className="gap-1 justify-start"
                  onClick={() => {
                    setSelectedTemplate(template.id);
                    setGenerationPrompt(template.description);
                  }}
                >
                  {template.icon}
                  <span className="truncate">{template.name}</span>
                </Button>
              ))}
            </div>

            <div className="flex gap-4 items-center">
              <div className="flex items-center gap-2">
                <Switch
                  checked={biologicalMode}
                  onCheckedChange={setBiologicalMode}
                />
                <label className="text-sm">Biological Patterns</label>
              </div>
              <div className="flex items-center gap-2">
                <Switch
                  checked={quantumMode}
                  onCheckedChange={setQuantumMode}
                />
                <label className="text-sm">Quantum Ready</label>
              </div>
              <div className="flex items-center gap-2">
                <Switch
                  checked={ethicsEnforced}
                  onCheckedChange={setEthicsEnforced}
                />
                <label className="text-sm">Ethics Enforcement</label>
              </div>
            </div>
          </div>
        </CardContent>
      </Card>

      {/* Main IDE */}
      <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
        <div className="lg:col-span-2">
          <Tabs value={activeTab} onValueChange={setActiveTab}>
            <TabsList className="grid w-full grid-cols-3">
              <TabsTrigger value="editor">Code Editor</TabsTrigger>
              <TabsTrigger value="visual">Visual Builder</TabsTrigger>
              <TabsTrigger value="output">Output</TabsTrigger>
            </TabsList>

            <TabsContent value="editor" className="space-y-4">
              <Card className="rounded-2xl shadow-sm">
                <CardHeader>
                  <div className="flex items-center justify-between">
                    <CardTitle className="text-base flex items-center gap-2">
                      <FileCode className="h-4 w-4" />
                      Algorithm Code
                    </CardTitle>
                    <div className="flex gap-2">
                      <Select value={language} onValueChange={(v: any) => setLanguage(v)}>
                        <SelectTrigger className="w-32">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="python">Python</SelectItem>
                          <SelectItem value="typescript">TypeScript</SelectItem>
                          <SelectItem value="rust">Rust</SelectItem>
                        </SelectContent>
                      </Select>
                      <Button size="sm" variant="outline" className="gap-1">
                        <Copy className="h-3 w-3" /> Copy
                      </Button>
                      <Button size="sm" variant="outline" className="gap-1">
                        <Save className="h-3 w-3" /> Save
                      </Button>
                      <Button 
                        size="sm" 
                        variant="default" 
                        className="gap-1"
                        onClick={runAlgorithm}
                        disabled={isRunning}
                      >
                        {isRunning ? (
                          <>
                            <Pause className="h-3 w-3" /> Stop
                          </>
                        ) : (
                          <>
                            <Play className="h-3 w-3" /> Run
                          </>
                        )}
                      </Button>
                    </div>
                  </div>
                </CardHeader>
                <CardContent>
                  <div className="h-[500px] border rounded-xl overflow-hidden">
                    <Editor
                      height="100%"
                      language={language}
                      value={code || "# Click 'Generate' to create an algorithm or select a template"}
                      theme="vs-dark"
                      options={{
                        minimap: { enabled: false },
                        fontSize: 14,
                        lineNumbers: "on",
                        scrollBeyondLastLine: false,
                        automaticLayout: true,
                        tabSize: 2,
                      }}
                      onChange={(value) => setCode(value || "")}
                    />
                  </div>
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="visual">
              <VisualAlgorithmBuilder onGenerate={setCode} />
            </TabsContent>

            <TabsContent value="output">
              <Card className="rounded-2xl shadow-sm">
                <CardHeader>
                  <CardTitle className="text-base flex items-center gap-2">
                    <Terminal className="h-4 w-4" />
                    Execution Output
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="h-[500px] bg-black text-green-400 font-mono text-sm p-4 rounded-xl overflow-auto">
                    <pre>{output || "No output yet. Run an algorithm to see results."}</pre>
                  </div>
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>
        </div>

        {/* Side Panels */}
        <div className="space-y-4">
          <PerformanceMonitor metrics={performanceMetrics} />

          <Card className="rounded-2xl shadow-sm">
            <CardHeader>
              <CardTitle className="text-base flex items-center gap-2">
                <Brain className="h-4 w-4" />
                AI Assistant
              </CardTitle>
            </CardHeader>
            <CardContent className="space-y-3">
              <Alert>
                <Lightbulb className="h-4 w-4" />
                <AlertDescription>
                  Your algorithm uses the Crowe Framework for multi-objective optimization. Consider adding parallel processing for better performance.
                </AlertDescription>
              </Alert>
              <div className="space-y-2">
                <Button variant="outline" size="sm" className="w-full justify-start gap-2">
                  <Zap className="h-3 w-3" />
                  Optimize for Speed
                </Button>
                <Button variant="outline" size="sm" className="w-full justify-start gap-2">
                  <Target className="h-3 w-3" />
                  Improve Accuracy
                </Button>
                <Button variant="outline" size="sm" className="w-full justify-start gap-2">
                  <Shield className="h-3 w-3" />
                  Add Safety Checks
                </Button>
                <Button variant="outline" size="sm" className="w-full justify-start gap-2">
                  <Database className="h-3 w-3" />
                  Connect to Database
                </Button>
              </div>
            </CardContent>
          </Card>

          <Card className="rounded-2xl shadow-sm">
            <CardHeader>
              <CardTitle className="text-base flex items-center gap-2">
                <Package className="h-4 w-4" />
                Algorithm Library
              </CardTitle>
            </CardHeader>
            <CardContent>
              <ScrollArea className="h-64">
                <div className="space-y-2">
                  {["Tanimoto Similarity", "MACCS Keys", "Morgan Fingerprints", "Pharmacophore Matching", "Shape Overlay", "Electrostatic Similarity"].map(algo => (
                    <Button
                      key={algo}
                      variant="ghost"
                      size="sm"
                      className="w-full justify-start text-xs"
                      onClick={() => setGenerationPrompt(`Implement ${algo} algorithm`)}
                    >
                      <Fingerprint className="h-3 w-3 mr-2" />
                      {algo}
                    </Button>
                  ))}
                </div>
              </ScrollArea>
            </CardContent>
          </Card>

          <Card className="rounded-2xl shadow-sm">
            <CardHeader>
              <CardTitle className="text-base flex items-center gap-2">
                <Bot className="h-4 w-4" />
                Active Agents
              </CardTitle>
            </CardHeader>
            <CardContent className="space-y-2">
              <div className="flex items-center justify-between">
                <div className="flex items-center gap-2">
                  <div className="h-2 w-2 rounded-full bg-green-500 animate-pulse" />
                  <span className="text-sm">Algorithm Generator</span>
                </div>
                <Badge variant="outline" className="text-xs">Active</Badge>
              </div>
              <div className="flex items-center justify-between">
                <div className="flex items-center gap-2">
                  <div className="h-2 w-2 rounded-full bg-blue-500 animate-pulse" />
                  <span className="text-sm">Code Optimizer</span>
                </div>
                <Badge variant="outline" className="text-xs">Ready</Badge>
              </div>
              <div className="flex items-center justify-between">
                <div className="flex items-center gap-2">
                  <div className="h-2 w-2 rounded-full bg-purple-500" />
                  <span className="text-sm">Ethics Validator</span>
                </div>
                <Badge variant="outline" className="text-xs">Idle</Badge>
              </div>
            </CardContent>
          </Card>
        </div>
      </div>

      {/* Footer */}
      <motion.div 
        initial={{ opacity: 0, y: 8 }} 
        animate={{ opacity: 1, y: 0 }} 
        transition={{ duration: 0.4 }} 
        className="text-center text-xs text-muted-foreground pt-6"
      >
        CriOS Algorithm Studio • AI-Powered Drug Discovery • © 2025 Crowe BioSystems
      </motion.div>
    </div>
  );
}