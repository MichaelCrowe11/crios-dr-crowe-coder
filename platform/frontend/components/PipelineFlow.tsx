// Pipeline Flow Visualization Component
// Shows drug discovery pipeline stages and progress

import React from 'react';
import { 
  Box, 
  Paper, 
  Typography, 
  Stepper, 
  Step, 
  StepLabel, 
  StepContent,
  Chip,
  LinearProgress,
  Avatar,
  AvatarGroup
} from '@mui/material';
import {
  CheckCircle as CheckIcon,
  RadioButtonUnchecked as PendingIcon,
  Loop as ProcessingIcon
} from '@mui/icons-material';
import { motion } from 'framer-motion';

interface PipelineStage {
  name: string;
  task: string;
  agents: string[];
  status: 'pending' | 'processing' | 'completed';
  progress: number;
  compounds: number;
}

interface PipelineFlowProps {
  pipeline: string;
  stages: PipelineStage[];
  currentStage: number;
}

const PipelineFlow: React.FC<PipelineFlowProps> = ({ pipeline, stages, currentStage }) => {
  const getStageIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CheckIcon color="success" />;
      case 'processing':
        return <ProcessingIcon color="primary" className="animate-spin" />;
      default:
        return <PendingIcon color="disabled" />;
    }
  };

  const getStageColor = (status: string) => {
    switch (status) {
      case 'completed':
        return 'success';
      case 'processing':
        return 'primary';
      default:
        return 'default';
    }
  };

  return (
    <Box>
      <Typography variant="h5" gutterBottom>
        {pipeline.replace('_', ' ').toUpperCase()} Pipeline
      </Typography>

      <Stepper activeStep={currentStage} orientation="vertical">
        {stages.map((stage, index) => (
          <Step key={index} completed={stage.status === 'completed'}>
            <StepLabel
              StepIconComponent={() => getStageIcon(stage.status)}
              optional={
                <Typography variant="caption">
                  {stage.agents.length} agents • {stage.compounds} compounds
                </Typography>
              }
            >
              <Typography variant="h6">{stage.name}</Typography>
              <Typography variant="body2" color="textSecondary">
                {stage.task.replace('_', ' ').charAt(0).toUpperCase() + stage.task.slice(1).replace('_', ' ')}
              </Typography>
            </StepLabel>
            <StepContent>
              <motion.div
                initial={{ opacity: 0, y: -10 }}
                animate={{ opacity: 1, y: 0 }}
                transition={{ duration: 0.3 }}
              >
                <Paper sx={{ p: 2, mb: 2 }}>
                  {/* Progress Bar */}
                  {stage.status === 'processing' && (
                    <Box sx={{ mb: 2 }}>
                      <Box sx={{ display: 'flex', justifyContent: 'space-between', mb: 1 }}>
                        <Typography variant="body2">Progress</Typography>
                        <Typography variant="body2">{stage.progress}%</Typography>
                      </Box>
                      <LinearProgress 
                        variant="determinate" 
                        value={stage.progress} 
                        sx={{ height: 8, borderRadius: 4 }}
                      />
                    </Box>
                  )}

                  {/* Agent Avatars */}
                  <Box sx={{ mb: 2 }}>
                    <Typography variant="body2" gutterBottom>
                      Active Agents:
                    </Typography>
                    <AvatarGroup max={5} sx={{ justifyContent: 'flex-start' }}>
                      {stage.agents.map((agent, i) => (
                        <Avatar 
                          key={i}
                          sx={{ 
                            width: 32, 
                            height: 32,
                            bgcolor: `hsl(${i * 60}, 70%, 50%)`
                          }}
                        >
                          {agent[0]}
                        </Avatar>
                      ))}
                    </AvatarGroup>
                  </Box>

                  {/* Stage Metrics */}
                  <Box sx={{ display: 'flex', gap: 1, flexWrap: 'wrap' }}>
                    <Chip 
                      label={`${stage.compounds} compounds`} 
                      size="small" 
                      color={getStageColor(stage.status) as any}
                    />
                    {stage.status === 'processing' && (
                      <Chip 
                        label="Running" 
                        size="small" 
                        color="primary"
                        variant="outlined"
                      />
                    )}
                    {stage.status === 'completed' && (
                      <Chip 
                        label="Completed" 
                        size="small" 
                        color="success"
                        variant="outlined"
                      />
                    )}
                  </Box>
                </Paper>
              </motion.div>
            </StepContent>
          </Step>
        ))}
      </Stepper>

      {/* Pipeline Summary */}
      {currentStage === stages.length && (
        <motion.div
          initial={{ opacity: 0, scale: 0.9 }}
          animate={{ opacity: 1, scale: 1 }}
          transition={{ duration: 0.5 }}
        >
          <Paper sx={{ p: 3, mt: 3, bgcolor: 'success.dark' }}>
            <Typography variant="h6" gutterBottom>
              Pipeline Complete! 
            </Typography>
            <Typography variant="body1">
              Total compounds generated: {stages.reduce((sum, s) => sum + s.compounds, 0)}
            </Typography>
            <Typography variant="body2" color="textSecondary">
              Total agents deployed: {new Set(stages.flatMap(s => s.agents)).size}
            </Typography>
          </Paper>
        </motion.div>
      )}
    </Box>
  );
};

export default PipelineFlow;