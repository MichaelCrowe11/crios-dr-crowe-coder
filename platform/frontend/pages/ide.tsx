// CriOS Immersive IDE - Standalone Page
// Full-screen development environment for algorithms and compounds

import React from 'react';
import Head from 'next/head';
import { 
  ThemeProvider,
  createTheme,
  CssBaseline
} from '@mui/material';
import dynamic from 'next/dynamic';

// Dynamic import to avoid SSR issues
const ImmersiveIDE = dynamic(() => import('../components/ImmersiveIDE'), { ssr: false });

// Dark theme optimized for coding
const ideTheme = createTheme({
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
      default: '#0d1117',
      paper: '#161b22'
    },
    text: {
      primary: '#f0f6fc',
      secondary: '#8b949e'
    }
  },
  typography: {
    fontFamily: '"JetBrains Mono", "Fira Code", "Monaco", monospace',
    fontSize: 13,
    body1: {
      fontFamily: '"Inter", "Roboto", sans-serif'
    }
  },
  components: {
    MuiPaper: {
      styleOverrides: {
        root: {
          backgroundImage: 'linear-gradient(rgba(255, 255, 255, 0.02), rgba(255, 255, 255, 0.02))',
          backdropFilter: 'blur(10px)',
          borderRadius: 8,
          border: '1px solid rgba(240, 246, 252, 0.1)'
        }
      }
    },
    MuiButton: {
      styleOverrides: {
        root: {
          borderRadius: 6,
          textTransform: 'none',
          fontWeight: 500
        }
      }
    }
  }
});

export default function IDEPage() {
  return (
    <ThemeProvider theme={ideTheme}>
      <CssBaseline />
      <Head>
        <title>CriOS IDE - Dr. Crowe Coder Development Environment</title>
        <meta name="description" content="Immersive development environment for biological computing algorithms and compound discovery" />
        <link rel="preconnect" href="https://fonts.googleapis.com" />
        <link rel="preconnect" href="https://fonts.gstatic.com" crossOrigin="anonymous" />
        <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;500;600;700&family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet" />
      </Head>
      
      <ImmersiveIDE />
    </ThemeProvider>
  );
}