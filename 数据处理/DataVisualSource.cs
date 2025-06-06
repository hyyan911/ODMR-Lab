﻿using ODMR_Lab.IO操作;
using ODMR_Lab.基本控件;
using ODMR_Lab.实验部分.ODMR实验.参数;
using ODMR_Lab.实验部分.场效应器件测量;
using ODMR_Lab.实验部分.磁场调节;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ODMR_Lab.数据处理
{
    /// <summary>
    /// 数据可视化对象
    /// </summary>
    public class DataVisualSource
    {
        /// <summary>
        /// 参数，分为标题和内容 
        /// </summary>
        public Dictionary<string, string> Params { get; set; } = new Dictionary<string, string>();

        /// <summary>
        /// 一维图表数据
        /// </summary>
        public List<ChartData1D> ChartDataSource1D { get; set; } = new List<ChartData1D>();

        /// <summary>
        /// 一维图表数据
        /// </summary>
        public List<ChartData2D> ChartDataSource2D { get; set; } = new List<ChartData2D>();

        /// <summary>
        /// 从文件加载
        /// </summary>
        /// <returns></returns>
        public static DataVisualSource LoadFromFile(string path)
        {
            ExperimentFileTypes type = ExperimentObject<ExpParamBase, ConfigBase>.GetExpType(path);
            if (type == ExperimentFileTypes.None) throw new Exception("此文件为不支持预览的类型,无法打开。\n文件路径:" + path);
            if (type == ExperimentFileTypes.源表IV测量数据)
            {
                IVMeasureExpObject obj = new IVMeasureExpObject();
                obj.ReadFromFile(path);
                return obj.ToDataVisualSource();
            }
            if (type == ExperimentFileTypes.自定义数据)
            {
                UserCustomExpObject obj = new UserCustomExpObject();
                obj.ReadFromFile(path);
                return obj.ToDataVisualSource();
            }
            if (type == ExperimentFileTypes.温度监测数据)
            {
                TemperatureExpObject obj = new TemperatureExpObject();
                obj.ReadFromFile(path);
                return obj.ToDataVisualSource();
            }
            if (type == ExperimentFileTypes.ODMR实验)
            {
                SequenceFileExpObject obj = new SequenceFileExpObject();
                obj.ReadFromFile(path);
                return obj.ToDataVisualSource();
            }
            throw new Exception("此文件为不支持预览的类型,无法打开");
        }
    }
}
