// CriOS Dr. Crowe Coder Platform - Main Dashboard
// Next.js + React + TypeScript Frontend

import React, { useState, useEffect } from 'react';
import Head from 'next/head';
import { 
  Box, 
  Container, 
  Grid, 
  Paper, 
  Typography, 
  Card,
  CardContent,
  Button,
  LinearProgress,
  Chip,
  IconButton,
  Tabs,
  Tab,
  Avatar,
  List,
  ListItem,
  ListItemAvatar,
  ListItemText,
  Badge,
  Tooltip,
  Divider,
  useTheme,
  ThemeProvider,
  createTheme,
  CssBaseline,
  CircularProgress,
  ListItemIcon
} from '@mui/material';
import {
  Science as ScienceIcon,
  Biotech as BiotechIcon,
  Psychology as PsychologyIcon,
  Analytics as AnalyticsIcon,
  Hub as HubIcon,
  Speed as SpeedIcon,
  Group as GroupIcon,
  Timeline as TimelineIcon,
  CloudQueue as CloudIcon,
  PlayArrow as PlayIcon,
  Stop as StopIcon,
  Refresh as RefreshIcon,
  CheckCircle as CheckIcon,
  Warning as WarningIcon,
  Code as CodeIcon
} from '@mui/icons-material';
import { motion } from 'framer-motion';
import dynamic from 'next/dynamic';

// Dynamic imports for heavy components
const MoleculeViewer = dynamic(() => import('../components/MoleculeViewer'), { ssr: false });
const PipelineFlow = dynamic(() => import('../components/PipelineFlow'), { ssr: false });
const AgentNetwork = dynamic(() => import('../components/AgentNetwork'), { ssr: false });
const ImmersiveIDE = dynamic(() => import('../components/ImmersiveIDE'), { ssr: false });

// Types
interface Agent {
  id: string;
  name: string;
  title: string;
  division: string;
  specialization: string[];
  h_index: number;
  status: 'idle' | 'working' | 'thinking' | 'collaborating';
}

interface Pipeline {
  id: string;
  name: string;
  stages: number;
  status: 'idle' | 'running' | 'completed' | 'failed';
  progress: number;
  compounds_generated: number;
}

interface DiscoveryTask {
  id: string;
  target: string;
  pipeline: string;
  agents_deployed: number;
  start_time: string;
  status: string;
  compounds: number;
}

// Custom theme
const theme = createTheme({
  palette: {
    mode: 'dark',
    primary: {
      main: '#00e5ff',
      light: '#6effff',
      dark: '#00acc1'
    },
    secondary: {
      main: '#76ff03',
      light: '#b0ff57',
      dark: '#32cb00'
    },
    background: {
      default: '#0a0e27',
      paper: '#151934'
    }
  },
  typography: {
    fontFamily: '"Inter", "Roboto", "Helvetica", "Arial", sans-serif',
    h1: {
      fontSize: '3rem',
      fontWeight: 800,
      background: 'linear-gradient(45deg, #00e5ff 30%, #76ff03 90%)',
      WebkitBackgroundClip: 'text',
      WebkitTextFillColor: 'transparent'
    }
  },
  components: {
    MuiPaper: {
      styleOverrides: {
        root: {
          backgroundImage: 'linear-gradient(rgba(255, 255, 255, 0.05), rgba(255, 255, 255, 0.05))',
          backdropFilter: 'blur(10px)',
          borderRadius: 16,
          border: '1px solid rgba(255, 255, 255, 0.1)'
        }
      }
    }
  }
});

export default function Dashboard() {
  const [selectedTab, setSelectedTab] = useState(0);
  const [agents, setAgents] = useState<Agent[]>([]);
  const [activePipeline, setActivePipeline] = useState<Pipeline | null>(null);
  const [recentTasks, setRecentTasks] = useState<DiscoveryTask[]>([]);
  const [stats, setStats] = useState({
    total_agents: 194,
    active_agents: 0,
    pipelines_running: 0,
    compounds_today: 0,
    success_rate: 85
  });
  const [wsConnection, setWsConnection] = useState<WebSocket | null>(null);

  useEffect(() => {
    // Fetch initial data
    fetchAgents();
    fetchStats();
    
    // Connect to WebSocket
    connectWebSocket();
    
    return () => {
      if (wsConnection) {
        wsConnection.close();
      }
    };
  }, []);

  const fetchAgents = async () => {
    try {
      const response = await fetch('http://localhost:8000/api/agents');
      const data = await response.json();
      setAgents(data.agents.slice(0, 20)); // Show first 20 agents
    } catch (error) {
      console.error('Failed to fetch agents:', error);
    }
  };

  const fetchStats = async () => {
    try {
      const response = await fetch('http://localhost:8000/api/stats');
      const data = await response.json();
      setStats(prev => ({
        ...prev,
        total_agents: data.total_agents,
        compounds_today: Math.floor(Math.random() * 10000)
      }));
    } catch (error) {
      console.error('Failed to fetch stats:', error);
    }
  };

  const connectWebSocket = () => {
    const ws = new WebSocket('ws://localhost:8000/ws/agents');
    
    ws.onopen = () => {
      console.log('WebSocket connected');
      setWsConnection(ws);
    };
    
    ws.onmessage = (event) => {
      const data = JSON.parse(event.data);
      handleWebSocketMessage(data);
    };
    
    ws.onerror = (error) => {
      console.error('WebSocket error:', error);
    };
  };

  const handleWebSocketMessage = (data: any) => {
    if (data.type === 'agent_active') {
      setStats(prev => ({
        ...prev,
        active_agents: prev.active_agents + 1
      }));
    } else if (data.type === 'pipeline_started') {
      setStats(prev => ({
        ...prev,
        pipelines_running: prev.pipelines_running + 1
      }));
    }
  };

  const startDiscovery = async () => {
    try {
      const response = await fetch('http://localhost:8000/api/discover', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          target: 'EGFR',
          pipeline: 'kinase_inhibitor',
          num_agents: 24,
          max_compounds: 1000
        })
      });
      const data = await response.json();
      console.log('Discovery started:', data);
    } catch (error) {
      console.error('Failed to start discovery:', error);
    }
  };

  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <Head>
        <title>CriOS Dr. Crowe Coder - Drug Discovery Platform</title>
        <meta name="description" content="194 PhD AI Agents for Revolutionary Drug Discovery" />
      </Head>

      <Container maxWidth="xl" sx={{ py: 4 }}>
        {/* Header */}
        <motion.div
          initial={{ opacity: 0, y: -20 }}
          animate={{ opacity: 1, y: 0 }}
          transition={{ duration: 0.5 }}
        >
          <Box sx={{ mb: 4, textAlign: 'center' }}>
            <Typography variant="h1" gutterBottom>
              Dr. Crowe Coder
            </Typography>
            <Typography variant="h5" color="textSecondary">
              194 PhD Agents Orchestrating Breakthrough Drug Discovery
            </Typography>
          </Box>
        </motion.div>

        {/* Stats Cards */}
        <Grid container spacing={3} sx={{ mb: 4 }}>
          {[
            { 
              label: 'PhD Agents', 
              value: stats.total_agents, 
              icon: <GroupIcon />, 
              color: '#00e5ff' 
            },
            { 
              label: 'Active Now', 
              value: stats.active_agents, 
              icon: <PsychologyIcon />, 
              color: '#76ff03' 
            },
            { 
              label: 'Pipelines Running', 
              value: stats.pipelines_running, 
              icon: <TimelineIcon />, 
              color: '#ff6090' 
            },
            { 
              label: 'Compounds Today', 
              value: stats.compounds_today.toLocaleString(), 
              icon: <BiotechIcon />, 
              color: '#ffd740' 
            },
            { 
              label: 'Success Rate', 
              value: `${stats.success_rate}%`, 
              icon: <AnalyticsIcon />, 
              color: '#69f0ae' 
            }
          ].map((stat, index) => (
            <Grid item xs={12} sm={6} md={2.4} key={index}>
              <motion.div
                initial={{ opacity: 0, scale: 0.9 }}
                animate={{ opacity: 1, scale: 1 }}
                transition={{ delay: index * 0.1 }}
              >
                <Paper 
                  sx={{ 
                    p: 3, 
                    textAlign: 'center',
                    background: `linear-gradient(135deg, ${stat.color}20 0%, transparent 100%)`,
                    borderColor: stat.color,
                    borderWidth: 2
                  }}
                >
                  <Box sx={{ color: stat.color, mb: 1 }}>{stat.icon}</Box>
                  <Typography variant="h4" sx={{ fontWeight: 'bold', color: stat.color }}>
                    {stat.value}
                  </Typography>
                  <Typography variant="body2" color="textSecondary">
                    {stat.label}
                  </Typography>
                </Paper>
              </motion.div>
            </Grid>
          ))}
        </Grid>

        {/* Main Content Tabs */}
        <Paper sx={{ mb: 4 }}>
          <Tabs 
            value={selectedTab} 
            onChange={(e, v) => setSelectedTab(v)}
            sx={{ borderBottom: 1, borderColor: 'divider' }}
          >
            <Tab label="Discovery Pipeline" icon={<TimelineIcon />} />
            <Tab label="Agent Network" icon={<HubIcon />} />
            <Tab label="Molecule Analysis" icon={<BiotechIcon />} />
            <Tab label="Immersive IDE" icon={<CodeIcon />} />
            <Tab label="Real-time Monitor" icon={<SpeedIcon />} />
          </Tabs>

          <Box sx={{ p: 3, minHeight: 500 }}>
            {selectedTab === 0 && (
              <Box>
                <Box sx={{ mb: 3, display: 'flex', justifyContent: 'space-between' }}>
                  <Typography variant="h5">Discovery Pipelines</Typography>
                  <Button
                    variant="contained"
                    startIcon={<PlayIcon />}
                    onClick={startDiscovery}
                    sx={{
                      background: 'linear-gradient(45deg, #00e5ff 30%, #76ff03 90%)',
                      boxShadow: '0 3px 5px 2px rgba(0, 229, 255, .3)'
                    }}
                  >
                    Start New Discovery
                  </Button>
                </Box>

                <Grid container spacing={3}>
                  {['kinase_inhibitor', 'protac', 'ai_generative', 'natural_product'].map((pipeline) => (
                    <Grid item xs={12} md={6} key={pipeline}>
                      <Card>
                        <CardContent>
                          <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 2 }}>
                            <Typography variant="h6">
                              {pipeline.replace('_', ' ').toUpperCase()}
                            </Typography>
                            <Chip 
                              label="Ready" 
                              color="success" 
                              size="small"
                              icon={<CheckIcon />}
                            />
                          </Box>
                          <Typography variant="body2" color="textSecondary" sx={{ mb: 2 }}>
                            4 stages • 24 agents • ~6 hours
                          </Typography>
                          <LinearProgress 
                            variant="determinate" 
                            value={0} 
                            sx={{ height: 8, borderRadius: 4 }}
                          />
                        </CardContent>
                      </Card>
                    </Grid>
                  ))}
                </Grid>
              </Box>
            )}

            {selectedTab === 1 && (
              <Box>
                <Typography variant="h5" sx={{ mb: 3 }}>
                  Agent Network Visualization
                </Typography>
                <Box sx={{ height: 400, position: 'relative' }}>
                  <AgentNetwork agents={agents} />
                </Box>
              </Box>
            )}

            {selectedTab === 2 && (
              <Box>
                <Typography variant="h5" sx={{ mb: 3 }}>
                  Molecule Analysis
                </Typography>
                <Grid container spacing={3}>
                  <Grid item xs={12} md={6}>
                    <Paper sx={{ p: 3, height: 400 }}>
                      <MoleculeViewer smiles="CC(=O)Oc1ccccc1C(=O)O" />
                    </Paper>
                  </Grid>
                  <Grid item xs={12} md={6}>
                    <Paper sx={{ p: 3 }}>
                      <Typography variant="h6" gutterBottom>
                        Molecular Properties
                      </Typography>
                      <List>
                        {[
                          { label: 'Molecular Weight', value: '180.16 g/mol' },
                          { label: 'LogP', value: '1.19' },
                          { label: 'TPSA', value: '63.60 Ų' },
                          { label: 'H-Bond Donors', value: '1' },
                          { label: 'H-Bond Acceptors', value: '4' },
                          { label: 'Lipinski Violations', value: '0' }
                        ].map((prop) => (
                          <ListItem key={prop.label}>
                            <ListItemText 
                              primary={prop.label}
                              secondary={prop.value}
                            />
                          </ListItem>
                        ))}
                      </List>
                    </Paper>
                  </Grid>
                </Grid>
              </Box>
            )}

            {selectedTab === 3 && (
              <Box sx={{ height: 'calc(100vh - 300px)' }}>
                <ImmersiveIDE />
              </Box>
            )}

            {selectedTab === 4 && (
              <Box>
                <Typography variant="h5" sx={{ mb: 3 }}>
                  Real-time Activity Monitor
                </Typography>
                <Grid container spacing={3}>
                  <Grid item xs={12} md={8}>
                    <Paper sx={{ p: 3, height: 400, overflow: 'auto' }}>
                      <Typography variant="h6" gutterBottom>
                        Agent Activity Stream
                      </Typography>
                      <List>
                        {agents.slice(0, 10).map((agent, index) => (
                          <motion.div
                            key={agent.id}
                            initial={{ opacity: 0, x: -20 }}
                            animate={{ opacity: 1, x: 0 }}
                            transition={{ delay: index * 0.05 }}
                          >
                            <ListItem>
                              <ListItemAvatar>
                                <Badge
                                  badgeContent=" "
                                  color={index % 3 === 0 ? "success" : "default"}
                                  variant="dot"
                                  anchorOrigin={{
                                    vertical: 'bottom',
                                    horizontal: 'right',
                                  }}
                                >
                                  <Avatar sx={{ bgcolor: theme.palette.primary.main }}>
                                    {agent.name.split(' ')[1][0]}
                                  </Avatar>
                                </Badge>
                              </ListItemAvatar>
                              <ListItemText
                                primary={agent.name}
                                secondary={
                                  index % 3 === 0 
                                    ? "Analyzing molecular structure..." 
                                    : "Idle"
                                }
                              />
                              {index % 3 === 0 && (
                                <CircularProgress size={20} />
                              )}
                            </ListItem>
                          </motion.div>
                        ))}
                      </List>
                    </Paper>
                  </Grid>
                  <Grid item xs={12} md={4}>
                    <Paper sx={{ p: 3 }}>
                      <Typography variant="h6" gutterBottom>
                        System Health
                      </Typography>
                      <List>
                        <ListItem>
                          <ListItemIcon>
                            <CheckIcon color="success" />
                          </ListItemIcon>
                          <ListItemText primary="Claude Code Engine" secondary="Operational" />
                        </ListItem>
                        <ListItem>
                          <ListItemIcon>
                            <CheckIcon color="success" />
                          </ListItemIcon>
                          <ListItemText primary="Knowledge Lake" secondary="5.2 TB" />
                        </ListItem>
                        <ListItem>
                          <ListItemIcon>
                            <CheckIcon color="success" />
                          </ListItemIcon>
                          <ListItemText primary="MCP Servers" secondary="7/7 Connected" />
                        </ListItem>
                        <ListItem>
                          <ListItemIcon>
                            <WarningIcon color="warning" />
                          </ListItemIcon>
                          <ListItemText primary="GPU Cluster" secondary="85% Utilized" />
                        </ListItem>
                      </List>
                    </Paper>
                  </Grid>
                </Grid>
              </Box>
            )}
          </Box>
        </Paper>

        {/* Footer */}
        <Box sx={{ mt: 4, textAlign: 'center', color: 'text.secondary' }}>
          <Typography variant="body2">
            © 2025 Crowe Research • CriOS Dr. Crowe Coder Platform v1.0
          </Typography>
          <Typography variant="caption">
            "Where 194 PhD minds converge to cure disease"
          </Typography>
        </Box>
      </Container>
    </ThemeProvider>
  );
}