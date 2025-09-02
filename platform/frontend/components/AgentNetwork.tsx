// Agent Network Visualization Component
// Interactive D3.js network graph showing 194 PhD agents

import React, { useEffect, useRef } from 'react';
import * as d3 from 'd3';
import { Box, useTheme } from '@mui/material';

interface Agent {
  id: string;
  name: string;
  division: string;
  status: string;
}

interface AgentNetworkProps {
  agents: Agent[];
}

const AgentNetwork: React.FC<AgentNetworkProps> = ({ agents }) => {
  const svgRef = useRef<SVGSVGElement>(null);
  const theme = useTheme();

  useEffect(() => {
    if (!svgRef.current || agents.length === 0) return;

    const width = 800;
    const height = 400;

    // Clear previous content
    d3.select(svgRef.current).selectAll('*').remove();

    const svg = d3.select(svgRef.current)
      .attr('width', '100%')
      .attr('height', '100%')
      .attr('viewBox', `0 0 ${width} ${height}`);

    // Create force simulation
    const simulation = d3.forceSimulation(agents as any)
      .force('charge', d3.forceManyBody().strength(-50))
      .force('center', d3.forceCenter(width / 2, height / 2))
      .force('collision', d3.forceCollide().radius(30))
      .force('x', d3.forceX(width / 2).strength(0.1))
      .force('y', d3.forceY(height / 2).strength(0.1));

    // Division colors
    const divisionColors: Record<string, string> = {
      molecular_design: '#00e5ff',
      biological_systems: '#76ff03',
      clinical_research: '#ff6090',
      data_science: '#ffd740',
      synthesis_chemistry: '#69f0ae',
      pharmacology: '#b388ff',
      regulatory_affairs: '#ff8a65',
      innovation_strategy: '#81d4fa',
      crowe_logic_core: '#ffffff'
    };

    // Create gradient definitions
    const defs = svg.append('defs');
    
    Object.entries(divisionColors).forEach(([division, color]) => {
      const gradient = defs.append('radialGradient')
        .attr('id', `gradient-${division}`)
        .attr('cx', '30%')
        .attr('cy', '30%');
      
      gradient.append('stop')
        .attr('offset', '0%')
        .attr('stop-color', color)
        .attr('stop-opacity', 0.8);
      
      gradient.append('stop')
        .attr('offset', '100%')
        .attr('stop-color', color)
        .attr('stop-opacity', 0.2);
    });

    // Create links (connect agents in same division)
    const links: any[] = [];
    agents.forEach((agent1, i) => {
      agents.forEach((agent2, j) => {
        if (i < j && agent1.division === agent2.division) {
          links.push({ source: agent1, target: agent2 });
        }
      });
    });

    // Add links
    const link = svg.append('g')
      .selectAll('line')
      .data(links)
      .enter()
      .append('line')
      .attr('stroke', '#ffffff')
      .attr('stroke-opacity', 0.1)
      .attr('stroke-width', 1);

    // Add nodes
    const node = svg.append('g')
      .selectAll('g')
      .data(agents)
      .enter()
      .append('g')
      .attr('cursor', 'pointer');

    // Add circles
    node.append('circle')
      .attr('r', 20)
      .attr('fill', (d: any) => `url(#gradient-${d.division})`)
      .attr('stroke', (d: any) => divisionColors[d.division] || '#ffffff')
      .attr('stroke-width', 2)
      .attr('stroke-opacity', 0.8);

    // Add status indicator
    node.append('circle')
      .attr('r', 5)
      .attr('cx', 15)
      .attr('cy', -15)
      .attr('fill', (d: any) => {
        if (d.status === 'working') return '#76ff03';
        if (d.status === 'thinking') return '#ffd740';
        return '#666666';
      });

    // Add labels
    node.append('text')
      .text((d: any) => d.name.split(' ')[2]?.[0] || d.name[0])
      .attr('text-anchor', 'middle')
      .attr('dy', '.35em')
      .attr('fill', '#ffffff')
      .attr('font-size', '12px')
      .attr('font-weight', 'bold');

    // Add hover effect
    node.on('mouseover', function(event, d: any) {
      d3.select(this).select('circle')
        .transition()
        .duration(200)
        .attr('r', 25);
    })
    .on('mouseout', function(event, d: any) {
      d3.select(this).select('circle')
        .transition()
        .duration(200)
        .attr('r', 20);
    });

    // Update positions on tick
    simulation.on('tick', () => {
      link
        .attr('x1', (d: any) => d.source.x)
        .attr('y1', (d: any) => d.source.y)
        .attr('x2', (d: any) => d.target.x)
        .attr('y2', (d: any) => d.target.y);

      node.attr('transform', (d: any) => `translate(${d.x},${d.y})`);
    });

    // Drag functionality
    const drag = d3.drag()
      .on('start', (event, d: any) => {
        if (!event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
      })
      .on('drag', (event, d: any) => {
        d.fx = event.x;
        d.fy = event.y;
      })
      .on('end', (event, d: any) => {
        if (!event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
      });

    node.call(drag as any);

    return () => {
      simulation.stop();
    };
  }, [agents, theme]);

  return (
    <Box sx={{ width: '100%', height: '100%', position: 'relative' }}>
      <svg ref={svgRef} style={{ width: '100%', height: '100%' }} />
    </Box>
  );
};

export default AgentNetwork;